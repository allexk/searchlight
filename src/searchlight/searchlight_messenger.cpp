/*
 * Copyright 2014, Brown University, Providence, RI.
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
 * @file searchlight_messenger.cpp
 * The implementation of the Searchlight Messenger.
 *
 * @author Alexander Kalinin
 */

#include "searchlight_messenger.h"
#include "searchlight_messages.pb.h"

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.messenger"));

namespace {

// SL message types
enum SearchlightMessageType {
    mtSLChunkRequest = 50, // 50 is a relatively large number to avoid clashes
    mtSLChunk
};

class SearchlightMessageDesc : public scidb::MessageDesc,
    private boost::noncopyable
{
 public:
    // Empty message (mtNone type)
    SearchlightMessageDesc() : scidb::MessageDesc() {}

    // Empty message (mtNone type) with a binary attachment
    SearchlightMessageDesc(boost::shared_ptr<scidb::SharedBuffer> binary) :
        scidb::MessageDesc(binary) {}

    // Is the message valid?
    virtual bool validate() override {
        if (MessageDesc::validate()) {
            return true;
        }

        switch (getMessageType()) {
            case mtSLChunkRequest:
            case mtSLChunk:
                break;
            default:
                return false;
        }

        /*
         * Not exactly true, since the base validate() might have returned
         * false because of the incorrect network protocol version. However,
         * we have no way to check that.
         */
        return true;
    }

protected:
    // Creates an appropriate protobuf record type
    virtual MessagePtr createRecord(MessageID messageType) override {
        switch (messageType) {
            case mtSLChunkRequest:
                return MessagePtr(new FetchChunk);
            case mtSLChunk:
                return MessagePtr(new scidb_msg::Chunk);
            default:
                LOG4CXX_ERROR(logger, "Unknown type of Searchlight message"
                        "to create!");
                return MessagePtr();
        }
    }
};

/**
 * This class is basically a wrapper around CompressedBuffer that is used for
 * holding data the buffer does not own.
 *
 * Why? We create such a buffer from the message to de-compress it in a chunk.
 * Unfortunately, CompressedBuffer frees its data, which is unacceptable in
 * our case -- the data still belongs to the message. One solution would be
 * to create a copy of the binary data and give it to the CompressedBuffer,
 * although it would be inefficient for large chunks.
 *
 * We have to override getData() for that, since these clowns from SciDb's team
 * call _virtual_ function free() from the destructor. We don't override
 * other functions, since the usage is really controlled and nothing bad
 * will happen with the data_ pointer.
 */
class CompressedDataHolder : public CompressedBuffer {
public:
    CompressedDataHolder(void *compressedData, int compressionMethod,
            size_t compressedSize, size_t decompressedSize) :
                CompressedBuffer(nullptr, compressionMethod,
                        compressedSize, decompressedSize),
                data_(compressedData) {}

    CompressedDataHolder() : CompressedBuffer(), data_(nullptr) {}

    virtual void *getData() const override {
        return data_;
    }

private:
    // Pointer to data we manage instead of CompressedBuf
    void *data_;
};

/*
 *  Helper function to retrieve CompressedBuffer from the MessageDesc.
 *  Why? Because DefaultMessageDescription does not contain an access method
 *  for the buffer and spits out an asio::const_buffer instead (yuck!).
 */
boost::shared_ptr<CompressedBuffer> RetrieveCompressedBuffer(
        const boost::shared_ptr<MessageDescription> &msg_desc,
        const boost::shared_ptr<scidb_msg::Chunk> &chunk_record) {
    // the buffer tinkering
    const boost::asio::const_buffer &buf = msg_desc->getBinary();
    const size_t compressed_size = boost::asio::buffer_size(buf);
    const void *compressed_data = boost::asio::buffer_cast<const void *>(buf);
    // the rest of the info
    const size_t decompressed_size = chunk_record->decompressed_size();
    const int compression_method = chunk_record->compression_method();

    if (!compressed_data || compressed_size == 0) {
        return boost::shared_ptr<CompressedDataHolder>();
    } else {
        // const_cast is safe here -- we're going to read and de-compress
        return boost::make_shared<CompressedDataHolder>(
                const_cast<void *>(compressed_data), compression_method,
                compressed_size, decompressed_size);
    }
}

}

void SearchlightMessenger::init(const boost::shared_ptr<Query> &query) {
    if (!initialized_) {
        query_ = query;

        // Establish query finalizer
        Query::Finalizer f =
                boost::bind(&SearchlightMessenger::deactivate, this, _1);
        query->pushFinalizer(f);

        // Establish message listeners
        RegisterHandlers(true);

        initialized_ = true;
    }
}

void SearchlightMessenger::deactivate(const boost::shared_ptr<Query> &) {
    if (initialized_) {
        reg_arrays_.clear();
        query_.reset();
        initialized_ = false;

        /*
         * We do not remove message handlers here. First of all, SciDb does not
         * have API for that -- handlers remain until termination. Secondly,
         * we just throw an exception if a request came out of context (query).
         *
         * An alternative would be to establish empty handlers, so SciDb
         * would drop messages by itself. However, we won't detect errors
         * that way.
         */
    }
}

void SearchlightMessenger::RegisterArray(const ArrayPtr &array) {
    reg_arrays_[array->getName()] = array;
}

uint32_t SearchlightMessenger::GetMessageSlot() {
    std::lock_guard<std::mutex> lock(mtx_);
    for (size_t i = 0; i < slots_.size(); i++) {
        if (!slots_[i].used_) {
            slots_[i].used_ = true;
            return i;
        }
    }

    // create a new slot
    slots_.emplace_back();
    return slots_.size() - 1;
}

void SearchlightMessenger::ReturnMessageSlot(size_t slot) {
    std::lock_guard<std::mutex> lock(mtx_);
    assert(slots_[slot].used_);
    slots_[slot].used_ = false;
    slots_[slot].data_ready_ = false;
    slots_[slot].data_ = nullptr;
}

void SearchlightMessenger::RegisterHandlers(bool init) {
    boost::shared_ptr<scidb::NetworkMessageFactory> factory =
            scidb::getNetworkMessageFactory();
    factory->addMessageType(mtSLChunkRequest,
            boost::bind(&SearchlightMessenger::CreateSLChunkRequest, this, _1),
            boost::bind(&SearchlightMessenger::HandleSLChunkRequest, this, _1));
    factory->addMessageType(mtSLChunk,
            boost::bind(&SearchlightMessenger::CreateSLChunk, this, _1),
            boost::bind(&SearchlightMessenger::HandleSLChunk, this, _1));
}

bool SearchlightMessenger::RequestChunk(InstanceID inst,
        const std::string &array_name, const Coordinates &pos, AttributeID attr,
        Chunk *chunk) {
    // check if the array is registered
    const auto &array_iter = reg_arrays_.find(array_name);
    if (array_iter == reg_arrays_.end()) {
        LOG4CXX_ERROR(logger, "Requesting chunk from unregistered array!");
        return false;
    }

    // message params
    const bool confirm_only = (chunk == nullptr);
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    const size_t slot_num = GetMessageSlot();
    MessageSlot &msg_slot = slots_[slot_num];

    std::shared_ptr<ChunkRequestData> chunk_req =
            std::make_shared<ChunkRequestData>();
    chunk_req->chunk_ = chunk;
    chunk_req->array_ = array_iter->second;
    msg_slot.data_ = chunk_req.get();
    msg_slot.data_ready_ = false;

    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg(new SearchlightMessageDesc);
    msg->initRecord(mtSLChunkRequest);
    msg->setQueryID(query->getQueryID());
    boost::shared_ptr<FetchChunk> record = msg->getRecord<FetchChunk>();

    // fill the record
    record->set_array_name(array_name);
    record->set_attribute_id(attr);
    record->set_confirm_only(confirm_only);
    for (size_t i = 0; i < pos.size(); i++) {
        record->add_position(pos[i]);
    }
    record->set_slot(slot_num);

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(inst, msg);

    // now we have to wait
    {
        std::unique_lock<std::mutex> slot_lock(msg_slot.mtx_);
        while (!msg_slot.data_ready_) {
            msg_slot.cond_.wait(slot_lock);
        }
        // safe to unlock since we write only once to the slot
    }

    // got the data
    const bool chunk_empty = chunk_req->empty_;
    ReturnMessageSlot(slot_num);

    return !chunk_empty;
}

MessagePtr SearchlightMessenger::CreateSLChunkRequest(MessageID id) {
    // cannot use the scidb Fetch message (it doesn't have coordinates field)
    return MessagePtr(new FetchChunk);
}

void SearchlightMessenger::HandleSLChunkRequest(
        const boost::shared_ptr<MessageDescription> &msg_desc) {

    // some correctness checking
    const scidb::QueryID query_id = msg_desc->getQueryId();
    const boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query_id != query->getQueryID()) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Query IDs on searchlight messengers do not match!";
    }

    // get the record
    boost::shared_ptr<FetchChunk> record =
            boost::dynamic_pointer_cast<FetchChunk>(msg_desc->getRecord());
    if (!record) {
        throw SYSTEM_EXCEPTION(scidb::SCIDB_SE_NETWORK,
                scidb::SCIDB_LE_INVALID_MESSAGE_FORMAT) <<
                        msg_desc->getMessageType();
    }

    // parse the record and fetch
    const std::string &array_name = record->array_name();
    const auto &reg_iter = reg_arrays_.find(array_name);
    if (reg_iter == reg_arrays_.end()) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Array is not registered on the remote messenger!";
    }

    const ArrayPtr &array = reg_iter->second;
    const AttributeID attr = record->attribute_id();
    const bool confirm_only = record->confirm_only();
    Coordinates pos(record->position_size());
    for (size_t i = 0; i < pos.size(); i++) {
        pos[i] = record->position(i);
    }

    // find the chunk
    // TODO: cache iterators
    auto array_iter = array->getConstIterator(attr);
    bool exists = array_iter->setPosition(pos);

    // create the reply record
    const bool attach_binary = exists && !confirm_only;
    boost::shared_ptr<scidb::CompressedBuffer> buffer;
    boost::shared_ptr<scidb::MessageDesc> msg;
    if (attach_binary) {
        buffer = boost::make_shared<scidb::CompressedBuffer>();
        msg = boost::make_shared<SearchlightMessageDesc>(buffer);
    } else {
        msg = boost::make_shared<SearchlightMessageDesc>();
    }

    msg->initRecord(mtSLChunk);
    msg->setQueryID(query->getQueryID());
    boost::shared_ptr<scidb_msg::Chunk> reply =
            msg->getRecord<scidb_msg::Chunk>();
    reply->set_obj_type(record->slot());
    reply->set_eof(!exists);
    if (attach_binary) {
        const ConstChunk &chunk = array_iter->getChunk();
        boost::shared_ptr<ConstRLEEmptyBitmap> empty_bitmap;
        if (chunk.isRLE() &&
                array->getArrayDesc().getEmptyBitmapAttribute() &&
                !chunk.getAttributeDesc().isEmptyIndicator()) {
            empty_bitmap = chunk.getEmptyBitmap();
        }
        chunk.compress(*buffer, empty_bitmap);
        empty_bitmap.reset(); // apparently, reset() is mandatory here (bug?)

        reply->set_sparse(chunk.isSparse());
        reply->set_rle(chunk.isRLE());
        reply->set_compression_method(buffer->getCompressionMethod());
        reply->set_decompressed_size(buffer->getDecompressedSize());
        reply->set_count(chunk.isCountKnown() ? chunk.count() : 0);
    }

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->sendMessage(msg_desc->getSourceInstanceID(), msg);
}

MessagePtr SearchlightMessenger::CreateSLChunk(MessageID id) {
    // we will re-use the system chunk message
    return MessagePtr(new scidb_msg::Chunk);
}

void SearchlightMessenger::HandleSLChunk(
        const boost::shared_ptr<MessageDescription> &msg_desc) {

    // some correctness checking
    const scidb::QueryID query_id = msg_desc->getQueryId();
    boost::shared_ptr<Query> query = Query::getValidQueryPtr(query_);
    if (query_id != query->getQueryID()) {
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Query IDs on searchlight messengers do not match!";
    }

    // get the record
    boost::shared_ptr<scidb_msg::Chunk> record =
          boost::dynamic_pointer_cast<scidb_msg::Chunk>(msg_desc->getRecord());
    if (!record) {
        throw SYSTEM_EXCEPTION(scidb::SCIDB_SE_NETWORK,
                scidb::SCIDB_LE_INVALID_MESSAGE_FORMAT) <<
                        msg_desc->getMessageType();
    }

    // determine slot
    const size_t slot_num = record->obj_type();
    assert(slot_num < slots_.size());
    MessageSlot &slot = slots_[slot_num];
    assert(slot.used_);
    ChunkRequestData *chunk_req = static_cast<ChunkRequestData *>(slot.data_);

    // lock and go!
    std::lock_guard<std::mutex> lock(slot.mtx_);
    chunk_req->empty_ = record->eof();

    // binary data, if any
    boost::shared_ptr<CompressedBuffer> compressed_buffer =
            RetrieveCompressedBuffer(msg_desc, record);

    if (compressed_buffer && chunk_req->chunk_) {
        Chunk *chunk = chunk_req->chunk_;

        chunk->setSparse(record->sparse());
        chunk->setRLE(record->rle());
        chunk->decompress(*compressed_buffer);
        chunk->setCount(record->count());
        assert(checkChunkMagic(*chunk));
        /*
         *  Note, we do not call chunk->write() here, since the caller might
         *  want to do something else. write() might cause unpin() for some
         *  types of chunks (e.g., LRU Chunk) and a potential swap-out. This
         *  would be inefficient. So, write()/unpin() is the caller's
         *  responsibility.
         */
    }

    slot.data_ready_ = true;
    slot.cond_.notify_one();
}

}
