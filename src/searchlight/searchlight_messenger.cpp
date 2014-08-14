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

#include <query/Operator.h>

namespace searchlight {

// The logger
static log4cxx::LoggerPtr logger(
        log4cxx::Logger::getLogger("searchlight.messenger"));

namespace {

// Creates a record of the approppriate type for the message
MessagePtr CreateMessageRecord(MessageID id) {
    switch (id) {
        case SearchlightMessenger::mtSLChunkRequest:
            return MessagePtr(new FetchChunk);
        case SearchlightMessenger::mtSLChunk:
            return MessagePtr(new scidb_msg::Chunk);
        case SearchlightMessenger::mtSLSolution:
            return MessagePtr(new SearchlightSolution);
        default:
            LOG4CXX_ERROR(logger, "Unknown type of Searchlight message"
                    "to create!");
            return MessagePtr();
    }
}

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
            case SearchlightMessenger::mtSLChunkRequest:
            case SearchlightMessenger::mtSLChunk:
            case SearchlightMessenger::mtSLSolution:
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
        return CreateMessageRecord(messageType);
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

std::ostream &operator<<(std::ostream &os, const Coordinates &coords) {
    os << '(';
    for (size_t i = 0; i < coords.size(); i++) {
        if (i > 0) {
            os << ", ";
        }
        os << coords[i];
    }
    os << ')';
    return os;
}

} /* namespace (unanimous) */

void SearchlightMessenger::RegisterQuery(
        const boost::shared_ptr<Query> &query) {
    std::lock_guard<std::mutex> lock(mtx_);

    const QueryID query_id = query->getQueryID();
    const auto &iter = served_queries_.find(query_id);
    if (iter == served_queries_.end()) {
        served_queries_.emplace(query_id,
                std::make_shared<QueryContext>(query));

        // Establish query finalizer
        Query::Finalizer f =
                boost::bind(&SearchlightMessenger::deactivate, this, _1);
        query->pushFinalizer(f);

        LOG4CXX_DEBUG(logger, "Registered query, id=" << query_id);
    }
}

void SearchlightMessenger::deactivate(const boost::shared_ptr<Query> &query) {
    std::lock_guard<std::mutex> lock(mtx_);

    const QueryID query_id = query->getQueryID();
    const auto &iter = served_queries_.find(query_id);
    assert(iter != served_queries_.end());

    std::ostringstream os;
    os << "Query statistics follows (id=" << query_id << "):\n";
    iter->second->stats_.print(os);
    logger->info(os.str(), LOG4CXX_LOCATION);

    served_queries_.erase(iter);
    LOG4CXX_DEBUG(logger, "Removed query, id=" << query_id);

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

SearchlightMessenger::QueryContextPtr SearchlightMessenger::GetQueryContext(
        QueryID query_id, bool lock) const {
    std::unique_lock<std::mutex> lockg;
    if (lock) {
        lockg = std::unique_lock<std::mutex>(mtx_);
    }

    const auto &iter = served_queries_.find(query_id);
    if (iter == served_queries_.end()) {
        LOG4CXX_ERROR(logger,
                "Query is not registered at messenger: id=" << query_id);
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Invalid messenger request: no query with this id";
    }
    assert(iter->second);
    return iter->second;
}

void SearchlightMessenger::RegisterArray(const boost::shared_ptr<Query> &query,
        const ArrayPtr &array) {
    QueryContextPtr query_ctx = GetQueryContext(query->getQueryID(), true);
    const std::string &array_name = array->getName();
    query_ctx->reg_arrays_.emplace(array_name,
            std::make_shared<QueryContext::ServedArray>(array));
    LOG4CXX_DEBUG(logger, "Registered array, qid=" << query->getQueryID()
            << ", array=" << array_name);
}

SearchlightMessenger::QueryContext::ServedArrayPtr
SearchlightMessenger::GetRegisteredArray(
        const SearchlightMessenger::QueryContextPtr &query_ctx,
        const std::string &name) const {
    const auto &iter = query_ctx->reg_arrays_.find(name);
    if (iter == query_ctx->reg_arrays_.end()) {
        LOG4CXX_ERROR(logger,
                "Requesting chunk from unregistered array: " << name);
        throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                << "Array is not registered on the messenger!";
    }
    return iter->second;
}

void SearchlightMessenger::RegisterUserMessageHandler(
        const boost::shared_ptr<Query> &query,
        SearchlightMessageType msg_type,
        const UserMessageHandler &handler) {
    QueryContextPtr query_ctx = GetQueryContext(query->getQueryID(), true);
    query_ctx->message_handlers_.emplace(msg_type, handler);
}

boost::shared_ptr<Query> SearchlightMessenger::GetRegisteredQuery(
        const SearchlightMessenger::QueryContextPtr &query_ctx) const {
    return Query::getValidQueryPtr(query_ctx->query_);
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
    using MessageCreator = scidb::NetworkMessageFactory::MessageCreator;

    boost::shared_ptr<scidb::NetworkMessageFactory> factory =
            scidb::getNetworkMessageFactory();
    factory->addMessageType(mtSLChunkRequest,
            MessageCreator(&CreateMessageRecord),
            boost::bind(&SearchlightMessenger::HandleSLChunkRequest, this, _1));
    factory->addMessageType(mtSLChunk,
            MessageCreator(&CreateMessageRecord),
            boost::bind(&SearchlightMessenger::HandleSLChunk, this, _1));

    // General messages
    factory->addMessageType(mtSLSolution,
            MessageCreator(&CreateMessageRecord),
            boost::bind(&SearchlightMessenger::HandleGeneralMessage, this, _1));
}

bool SearchlightMessenger::RequestChunk(const boost::shared_ptr<Query> &query,
        InstanceID inst, const std::string &array_name, const Coordinates &pos,
        AttributeID attr, Chunk *chunk) {
    // check if the array is registered
    QueryContextPtr query_ctx = GetQueryContext(query->getQueryID(), true);
    GetRegisteredQuery(query_ctx); // just error checking
    ArrayPtr array = GetRegisteredArray(query_ctx, array_name)->array_;

    // message params
    const bool confirm_only = (chunk == nullptr);
    const size_t slot_num = GetMessageSlot();
    MessageSlot &msg_slot = slots_[slot_num];

    std::shared_ptr<ChunkRequestData> chunk_req =
            std::make_shared<ChunkRequestData>();
    chunk_req->chunk_ = chunk;
    chunk_req->array_ = array;
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

    // log
    if (logger->isDebugEnabled()) {
        std::ostringstream os;
        os <<"Requesting chunk: qid=" << query->getQueryID() << ", name=" <<
                array_name << ", pos=" << pos << ", attr=" << attr <<
                ", data=" <<
                std::boolalpha << !confirm_only;
        logger->debug(os.str(), LOG4CXX_LOCATION);
    }

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(inst, msg);

    // now we have to wait
    {
        const auto wait_start_time = std::chrono::steady_clock::now();

        std::unique_lock<std::mutex> slot_lock(msg_slot.mtx_);
        query_ctx->stats_.msgs_sent_++;
        while (!msg_slot.data_ready_) {
            msg_slot.cond_.wait(slot_lock);
        }

        const auto wait_end_time = std::chrono::steady_clock::now();
        query_ctx->stats_.total_wait_time_ +=
                std::chrono::duration_cast<std::chrono::microseconds>(
                        wait_end_time - wait_start_time);

        // safe to unlock since we write only once to the slot
    }

    // got the data
    const bool chunk_empty = chunk_req->empty_;
    ReturnMessageSlot(slot_num);

    return !chunk_empty;
}

void SearchlightMessenger::HandleSLChunkRequest(
        const boost::shared_ptr<MessageDescription> &msg_desc) {

    // some correctness checking
    const QueryContextPtr query_ctx =
            GetQueryContext(msg_desc->getQueryId(), true);
    const boost::shared_ptr<Query> query = GetRegisteredQuery(query_ctx);

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
    const QueryContext::ServedArrayPtr serv_array =
            GetRegisteredArray(query_ctx, array_name);
    const ArrayPtr array = serv_array->array_;
    const AttributeID attr = record->attribute_id();
    const bool confirm_only = record->confirm_only();
    Coordinates pos(record->position_size());
    for (size_t i = 0; i < pos.size(); i++) {
        pos[i] = record->position(i);
    }

    // log
    if (logger->isDebugEnabled()) {
        std::ostringstream os;
        os <<"Got chunk request: qid=" << query->getQueryID() << ", name=" <<
                array_name << ", pos=" << pos << ", attr=" << attr <<
                ", data=" <<
                std::boolalpha << !confirm_only;
        logger->debug(os.str(), LOG4CXX_LOCATION);
    }

    // get iterator (and cache it)
    boost::shared_ptr<ConstArrayIterator> array_iter;
    {
        std::lock_guard<std::mutex> lock(mtx_);
        array_iter = serv_array->iters_[attr];
        if (!array_iter) {
            array_iter = array->getConstIterator(attr);
            serv_array->iters_[attr] = array_iter;
        }
        query_ctx->stats_.msgs_sent_++;
        query_ctx->stats_.msgs_received_++;
    }

    // find the chunk
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

        // stats
        {
            std::lock_guard<std::mutex> lock(mtx_);
            query_ctx->stats_.chunks_sent_++;
            query_ctx->stats_.chunk_data_sent_ += buffer->getSize();
        }
    }

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->sendMessage(msg_desc->getSourceInstanceID(), msg);
}

void SearchlightMessenger::HandleSLChunk(
        const boost::shared_ptr<MessageDescription> &msg_desc) {

    // some correctness checking
    const QueryContextPtr query_ctx =
            GetQueryContext(msg_desc->getQueryId(), true);
    const boost::shared_ptr<Query> query = GetRegisteredQuery(query_ctx);

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

    // Retrieve data
    chunk_req->empty_ = record->eof();

    // binary data, if any
    boost::shared_ptr<CompressedBuffer> compressed_buffer =
            RetrieveCompressedBuffer(msg_desc, record);

    if (compressed_buffer && chunk_req->chunk_) {
        // stats
        {
            std::lock_guard<std::mutex> lock(slot.mtx_);
            query_ctx->stats_.chunks_received_++;
            query_ctx->stats_.chunk_data_received_ +=
                    compressed_buffer->getSize();
        }

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

    // log
    if (logger->isDebugEnabled()) {
        std::ostringstream os;
        os <<"Got chunk: qid=" << query->getQueryID() << ", name=" <<
                chunk_req->array_->getName() << ", empty=" <<
                std::boolalpha << chunk_req->empty_;
        logger->debug(os.str(), LOG4CXX_LOCATION);
    }

    std::lock_guard<std::mutex> lock(slot.mtx_);
    query_ctx->stats_.msgs_received_++;
    slot.data_ready_ = true;
    slot.cond_.notify_one();
}

void SearchlightMessenger::HandleGeneralMessage(
        const boost::shared_ptr<MessageDescription> &msg_desc) {

    // get the query
    QueryContextPtr query_ctx;
    {
        std::lock_guard<std::mutex> lock(mtx_);
        query_ctx = GetQueryContext(msg_desc->getQueryId(), false);
        query_ctx->stats_.msgs_received_++;
    }
    const boost::shared_ptr<Query> query = GetRegisteredQuery(query_ctx);

    // get the record
    const MessagePtr record = msg_desc->getRecord();
    if (!record) {
        throw SYSTEM_EXCEPTION(scidb::SCIDB_SE_NETWORK,
                scidb::SCIDB_LE_INVALID_MESSAGE_FORMAT) <<
                        msg_desc->getMessageType();
    }

    const MessageID mid = msg_desc->getMessageType();
    const auto &iter = query_ctx->message_handlers_.find(mid);
    if (iter != query_ctx->message_handlers_.end()) {
        LOG4CXX_DEBUG(logger,
                "Delegating message to the user handler: id=" << mid);
        const InstanceID from_id = query->mapPhysicalToLogical(
                msg_desc->getSourceInstanceID());

        iter->second(from_id, record.get());
    } else {
        LOG4CXX_ERROR(logger,
                "Dropping message: cannot find handler: id=" << mid);
    }
}

void SearchlightMessenger::SendSolution(const boost::shared_ptr<Query> &query,
        bool eor,
        const std::vector<int64_t> &vals) const {

    // determine the coordinator
    const InstanceID coord_id = query->getCoordinatorID();
    if (coord_id == scidb::COORDINATOR_INSTANCE) {
        LOG4CXX_ERROR(logger, "Attempting to send a solution"
                "from the coordinator");
        return;
    }

    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg(new SearchlightMessageDesc);
    msg->initRecord(mtSLSolution);
    msg->setQueryID(query->getQueryID());
    boost::shared_ptr<SearchlightSolution> record =
            msg->getRecord<SearchlightSolution>();

    // fill the record
    record->set_eor(eor);
    if (!eor) {
        VarAssignment *asgn = record->mutable_solution();
        for (size_t i = 0; i < vals.size(); i++) {
            asgn->add_var_min(vals[i]);
            // We don't need var_max here
        }
    }

    // log
    if (logger->isDebugEnabled()) {
        std::ostringstream os;
        os <<"Sending solution to the coordinator: qid=" << query->getQueryID()
                << ", eor=" << std::boolalpha << eor;
        logger->debug(os.str(), LOG4CXX_LOCATION);
    }

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(coord_id, msg);

    // stats
    {
        std::lock_guard<std::mutex> lock(mtx_);
        const QueryContextPtr query_ctx =
                GetQueryContext(query->getQueryID(), false);
        query_ctx->stats_.msgs_sent_++;
    }
}

void SearchlightMessenger::Synchronize(const boost::shared_ptr<Query> &query) {
    LOG4CXX_DEBUG(logger, "Requesting SL Messenger synchronization...");
    // The same id should be okay here, even in case of subsequent barriers
    scidb::syncBarrier(0, query);
    LOG4CXX_DEBUG(logger, "Completed SL Messenger synchronization...");
}

}