/* Copyright 2014, Brown University, Providence, RI.
 *
 *                         All Rights Reserved
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose other than its incorporation into a
 * commercial product is hereby granted without fee, provided that the
 * above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of Brown University not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.
 *
 * BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
 * PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
 * ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/**
 * @file DepartArray.cpp
 *
 * Implementation of the DepartArray.
 *
 * @author Alexander Kalinin
 */

#include "DepartArray.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.depart_array"));

class DepartArrayIterator : public ConstArrayIterator {
public:
    DepartArrayIterator(const DepartArray &array, AttributeID attr) :
        array_(array),
        cache_(*array.cache_.get()),
        input_iter_(array.inputArray->getConstIterator(attr)),
        cache_iter_(array.cache_->array_.getIterator(attr)),
        attr_(attr),
        positioned_(false) {

        const scidb::Dimensions &dims = array.desc.getDimensions();
        for (const auto &dim: dims) {
            low_.push_back(dim.getStartMin()); // cannot be inf
            high_.push_back(dim.getEndMax());  // cannot be inf
            chunk_intervals_.push_back(dim.getChunkInterval());
        }
    }

    /**
     * Checks if the end of the array has been reached.
     *
     * @return true, if iterator has reached the end of the array
     */
    virtual bool end() override {
        if (!positioned_) {
            reset();
        }
        return current_[0] > high_[0];
    }

    /**
     * Position iterator to the next chunk in row-major order.
     */
    virtual void operator++() override {
        if (!positioned_) {
            reset();
        }
        if (end()) {
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_NO_CURRENT_CHUNK);
        }
        FindNextPosition();
    }

    /**
     * Returns current iterator position.
     *
     * @return current iterator position
     */
    virtual Coordinates const &getPosition() override {
        if (!positioned_) {
            reset();
        }
        if (end()) {
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_NO_CURRENT_CHUNK);
        }
        return current_;
    }

    /**
     * Sets iterator current position. If setPosition() fails, this iterator
     * becomes invalid and must be either reset() or positioned for further
     * use. end() will answer true for an invalid iterator.
     *
     * @return true if specified position is valid
     *     (belongs to the array and the chunk exists); false, otherwise
     */
    virtual bool setPosition(Coordinates const &pos) override {
        positioned_ = true;
        if (!array_.desc.contains(pos)) {
            Invalidate();
            return false;
        }

        Coordinates chunk_pos(pos);
        array_.desc.getChunkPositionFor(chunk_pos);
        current_ = chunk_pos;
        const bool exists = CheckPosition(current_);
        if (!exists) {
            Invalidate();
        }
        return exists;
    }

    /**
     * Resets iterator to the beginning.
     */
    virtual void reset() override {
        positioned_ = true;
        current_ = low_;
        if (!CheckPosition(current_)) {
            FindNextPosition();
        }
    }

    /**
     * Returns chunk at the current position.
     *
     * @return chunk at the current position
     */
    virtual const ConstChunk &getChunk() override {
        if (!positioned_) {
            reset();
        }
        if (end()) {
            throw USER_EXCEPTION(SCIDB_SE_EXECUTION, SCIDB_LE_NO_CURRENT_CHUNK);
        }
        return GetChunk();
    }

private:

    // Checks if the position is valid (i.e., the chunk at pos exists)
    bool CheckPosition(const Coordinates &pos);

    // Gets (and possibly loads data for) the chunk at the specified position
    // The chunk must exist at the current position
    const ConstChunk &GetChunk() const;

    // Checks if we need to set the bitmap for this chunk and sets it.
    // Sets the bitmap for the chunk and, optionally, fills the bitmap chunk.
    void CheckAndSetBitmapRLE(MemChunk *chunk) const;

    // Fetches a bitmap from the one of the remotes
    void FetchBitmapFromRemote(const Coordinates &pos) const;

    // Check if we already have the chunk in cache.
    bool CheckRemoteCache(const Address &addr) const {
        const auto &cached_chunk = cache_.remote_chunks_.find(addr);
        return cached_chunk != cache_.remote_chunks_.end();
    }

    // Increments current position
    void IncrementCurrent() {
        for (size_t i = current_.size() - 1; i >= 1; i--) {
            current_[i] += chunk_intervals_[i];
            if (current_[i] > high_[i]) {
                current_[i] = low_[i];
            } else {
                return;
            }
        }
        current_[0] += chunk_intervals_[0];
        // here we might fall out outside, which means end() == true
    }

    // Find next position
    void FindNextPosition() {
        IncrementCurrent();
        while (!end()) {
            if (CheckPosition(current_)) {
                break;
            }
            IncrementCurrent();
        }
    }

    // Makes iterator invalid. Means that end() will return true after that.
    void Invalidate() {
        current_ = high_;
        current_[0]++; // see end()
    }

    // The corresponding Depart Array
    const DepartArray &array_;

    // Cache of the array
    DepartArray::CacheMemArray &cache_;

    // Iterator over the input (no modifications)
    const boost::shared_ptr<ConstArrayIterator> input_iter_;

    // Iterator over the cache array (will modify the cache with remote chunks)
    const boost::shared_ptr<ArrayIterator> cache_iter_;

    // Low, high boundaries for the whole array and the current position
    std::vector<Coordinate> low_, high_, current_;

    // Chunk intervals (for faster ++ iteration)
    std::vector<int64_t> chunk_intervals_;

    // The attribute we're traversing
    const AttributeID attr_;

    // Is it positioned?
    bool positioned_;
};

DepartArray::DepartArray(const ArrayDesc &desc, const ArrayPtr &input,
        const ArrayDistribution &input_distr,
        const boost::shared_ptr<Query> &query) :
    DelegateArray(desc, input, false),
    messenger_(SearchlightMessenger::getInstance()),
    query_(query),
    input_distr_(input_distr) {

    /*
     * We need atomic cache check and creation of the array if it's not
     * in the cache. For this we need to create a fake query. If the
     * array already exists, we dismiss the fake query. While this might
     * create some overhead, destroying a fake query is not that
     * cumbersome and is done only once per DepartArray creation.
     */
    const std::string &array_name = desc.getName();
    const boost::shared_ptr<Query> fake_query =
            Query::createFakeQuery(
                    query->mapLogicalToPhysical(query->getCoordinatorID()),
                    query->mapLogicalToPhysical(query->getInstanceID()),
                    query->getCoordinatorLiveness());
    const auto emplaced = cache_cache_.emplace(array_name, desc, fake_query);
    if (!emplaced.second) {
        // The array already existed in the cache
        Query::destroyFakeQuery(fake_query.get());
        LOG4CXX_DEBUG(logger, "Array " << array_name <<
                " was found in the cache");
    }
    cache_ = emplaced.first;

    // Set input partitioning schema (at delegate array)
    this->desc.setPartitioningSchema(input_distr_.getPartitioningSchema());
}

boost::shared_ptr<ConstArrayIterator>
    DepartArray::getConstIterator(AttributeID attr) const {

    return boost::make_shared<DepartArrayIterator>(*this, attr);
}

bool DepartArrayIterator::CheckPosition(const Coordinates &pos) {
    boost::shared_ptr<Query> query = Query::getValidQueryPtr(array_.query_);
    const InstanceID local_instance = query->getInstanceID();

    // where is this chunk?
    InstanceID chunk_instance =
        scidb::getInstanceForChunk(query, pos, array_.desc,
                array_.input_distr_.getPartitioningSchema(),
                boost::shared_ptr<scidb::DistributionMapper>(),
                0, local_instance);

    if (chunk_instance == local_instance) {
        // Check the input array
        return input_iter_->setPosition(pos);
    } else {
        // The hard case: have to check remote chunks and maybe do a request
        std::unique_lock<std::mutex> lock(cache_.mtx_);
        const Address addr(attr_, pos);
        const bool chunk_fetched = CheckRemoteCache(addr);
        if (!chunk_fetched) {
            // Have to go to the remote; check if there's a request
            auto req = cache_.current_requests_.find(addr);
            bool repeat_request = true;
            while (repeat_request) {
                if (req != cache_.current_requests_.end()) {
                    // Already have a request -- wait
                    while (req != cache_.current_requests_.end()) {
                        cache_.cond_.wait(lock);
                        req = cache_.current_requests_.find(addr);
                    }
                    // Check if we need to repeat the request
                    repeat_request =
                            cache_.repeat_requests_.find(addr) !=
                                    cache_.repeat_requests_.end();
                } else {
                    // No request yet or repeat -- we're going to fetch
                    cache_.current_requests_.insert(addr);
                    cache_.repeat_requests_.erase(addr);
                    lock.unlock();

                    // check the chunk via the messenger (no chunk data!)
                    bool chunk_exists;
                    try {
                        chunk_exists = array_.messenger_->RequestChunk(
                            query, chunk_instance, array_.desc.getName(), pos,
                            attr_, nullptr);
                    } catch (...) {
                        LOG4CXX_DEBUG(logger, "Cannot fulfill the current "
                                "non-data chunk request. Another will repeat.");
                        cache_.repeat_requests_.insert(addr);
                        cache_.current_requests_.erase(addr);
                        cache_.cond_.notify_all();
                        throw;
                    }

                    // Working with structures again
                    lock.lock();
                    cache_.remote_chunks_.emplace(addr, chunk_exists);

                    // erase via addr, not req, since req might've been invalidated
                    cache_.current_requests_.erase(addr);
                    cache_.cond_.notify_all();
                    return chunk_exists;
                }
            }
        }

        // Must have the record in cache by now (mutex is still locked)
        assert(CheckRemoteCache(addr));
        return cache_.remote_chunks_[addr];
    }
}

const ConstChunk &DepartArrayIterator::GetChunk() const {
    boost::shared_ptr<Query> query = Query::getValidQueryPtr(array_.query_);
    const InstanceID local_instance = query->getInstanceID();

    // where is this chunk?
    InstanceID chunk_instance =
        scidb::getInstanceForChunk(query, current_, array_.desc,
                array_.input_distr_.getPartitioningSchema(),
                boost::shared_ptr<scidb::DistributionMapper>(),
                0, local_instance);

    if (chunk_instance == local_instance) {
        // Check the input array
        input_iter_->setPosition(current_);
        return input_iter_->getChunk();
    } else {
        // have to lock, since somebody might be changing the array or cache
        std::unique_lock<std::mutex> lock(cache_.mtx_);
        assert(CheckRemoteCache(Address(attr_, current_)));
        assert(cache_.remote_chunks_[Address(attr_, current_)]);

        const Address addr(attr_, current_);
        auto req = cache_.current_requests_.find(addr);
        /*
         * A chunk is fetched only when it's in the cache array already and
         * there are no current or repeat requests.
         */
        const bool chunk_fetched = cache_iter_->setPosition(current_) &&
                req == cache_.current_requests_.end() &&
                cache_.repeat_requests_.find(addr) ==
                        cache_.repeat_requests_.end();
        if (chunk_fetched) {
            // have already fetched the data
            return cache_iter_->getChunk();
        } else {
            // have to fetch from the remote
            while (true) {
                if (req != cache_.current_requests_.end()) {
                    // Already have a request -- wait
                    while (req != cache_.current_requests_.end()) {
                        cache_.cond_.wait(lock);
                        req = cache_.current_requests_.find(addr);
                    }
                    // Check if we need to repeat the request
                    const bool repeat_request =
                            cache_.repeat_requests_.find(addr) !=
                                    cache_.repeat_requests_.end();
                    if (!repeat_request) {
                        // chunk should be fetched by now
                        cache_iter_->setPosition(current_);
                        return cache_iter_->getChunk();
                    }
                } else {
                    // No request yet or repeat -- we're going to fetch
                    cache_.current_requests_.insert(addr);
                    cache_.repeat_requests_.erase(addr);

                    /*
                     * We might already have the chunk initialized here, if one
                     * of the previous requests failed and we're doing a repeat.
                     * newChunk() should still work for this situation, though.
                     *
                     * Conversion is safe here since we explicitly use MemArray.
                     */
                    MemChunk *new_chunk = dynamic_cast<MemChunk *>(
                                &(cache_iter_->newChunk(current_)));
                    assert(new_chunk);
                    const ConstChunk *bitmap_chunk = new_chunk->getBitmapChunk();
                    const bool need_bitmap =
                            bitmap_chunk && bitmap_chunk->getSize() == 0;
                    lock.unlock();

                    // get the chunk from messenger
                    try {
                        array_.messenger_->RequestChunk(
                            query, chunk_instance, array_.desc.getName(),
                            current_, attr_, new_chunk);
                    } catch (...) {
                        LOG4CXX_DEBUG(logger, "Cannot fulfill the current "
                                "chunk fetch. Another will repeat.");
                        cache_.repeat_requests_.insert(addr);
                        cache_.current_requests_.erase(addr);
                        cache_.cond_.notify_all();
                        throw;
                    }

                    /*
                     * Okay, here we want to retrieve the bitmap chunk if our
                     * new chunk is not RLE (rare). If it's not, it does not
                     * contain a bitmap attached.
                     */
                    if (!new_chunk->isRLE() && need_bitmap) {
                        try {
                            /*
                             *  Fetching bitmap might be interrupted as well.
                             *  While the corresponding repeat request will
                             *  be created, we still need to notify this
                             *  chunk's clients.
                             *
                             *  This, however, might result in refetching the
                             *  current chunk.
                             */
                            FetchBitmapFromRemote(current_);
                        } catch (...) {
                            LOG4CXX_DEBUG(logger, "Cannot fulfill the current "
                                  "chunk's bitmap fetch. Another will repeat.");
                            cache_.repeat_requests_.insert(addr);
                            cache_.current_requests_.erase(addr);
                            cache_.cond_.notify_all();
                            throw;
                        }
                    }

                    // Working with structures again
                    lock.lock();
                    if (new_chunk->isRLE() &&
                            !new_chunk->getAttributeDesc().isEmptyIndicator()) {
                        CheckAndSetBitmapRLE(new_chunk);
                    }
                    new_chunk->write(query); // this will un-pin it

                    // erase via addr, not req, since req might've been invalidated
                    cache_.current_requests_.erase(addr);
                    cache_.cond_.notify_all();
                    return *new_chunk;
                }
            }
        }
    }
}

void DepartArrayIterator::CheckAndSetBitmapRLE(MemChunk *chunk) const {
    assert(chunk->isRLE());
    // Removing const from the pointer is safe here
    ConstChunk * const bitmap_chunk =
            const_cast<ConstChunk *>(chunk->getBitmapChunk());
    if (bitmap_chunk && bitmap_chunk->getSize() == 0) {
        // Retrieve the bitmap; it must be here by design!
        const size_t bitmap_size = chunk->getBitmapSize();
        const size_t bitmap_offset = chunk->getSize() - bitmap_size;
        assert(bitmap_size > 0);

        const char * const bitmap_raw_data =
                static_cast<char *>(chunk->getData()) + bitmap_offset;
        const boost::shared_ptr<const ConstRLEEmptyBitmap> bitmap =
                make_shared<const ConstRLEEmptyBitmap>(bitmap_raw_data);

        /*
         *  Put it in the bitmap chunk. There is a potential for a data race
         *  here. Since somebody else might be doing bitmap request. In this
         *  case the request and this write will be served by different threads.
         *  However, these threads are going to write the same data. So, this
         *  should be safe.
         *
         *  FIXME: Check for the bitmap request in GetChunk() and either create
         *  the request for the bitmap chunk and fulfill it or wait until
         *  the "foreign" bitmap chunk request is fulfilled (since we need
         *  the bitmap for the data chunk to function properly).
         */
        Chunk *mod_bitmap_chunk =
                dynamic_cast<Chunk *>(bitmap_chunk);
        assert(mod_bitmap_chunk);
        mod_bitmap_chunk->pin();
        mod_bitmap_chunk->allocate(bitmap_size);
        bitmap->pack(static_cast<char *>(mod_bitmap_chunk->getDataForLoad()));
        mod_bitmap_chunk->setRLE(true);
        mod_bitmap_chunk->unPin();

        // Create fetch record
        cache_.remote_chunks_.emplace(
                Address(bitmap_chunk->getAttributeDesc().getId(),
                        bitmap_chunk->getFirstPosition(false)),
                true);
    }
}

void DepartArrayIterator::FetchBitmapFromRemote(const Coordinates &pos) const {
    const AttributeDesc *empty_attr = array_.desc.getEmptyBitmapAttribute();
    if (empty_attr) {
        auto iter = array_.getConstIterator(empty_attr->getId());
        if (iter->setPosition(pos)) {
            iter->getChunk();
        }
    }
}

void DepartArray::ClearPersistentCache() {
    LOG4CXX_INFO(logger, "Clearing cache for DepartArray...");
    cache_cache_.Clear();
}


// Cache of cache arrays (shared between all DepartArray clients)
DepartArray::CacheMemArrayCache DepartArray::cache_cache_;
}
