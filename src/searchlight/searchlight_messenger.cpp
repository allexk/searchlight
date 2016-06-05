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
        case SearchlightMessenger::mtSLControl:
            return MessagePtr(new SearchlightControl);
        case SearchlightMessenger::mtSLBalance:
            return MessagePtr(new SearchlightBalance);
        case SearchlightMessenger::mtSLMeta:
            return MessagePtr(new SearchlightMeta);
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
            case SearchlightMessenger::mtSLControl:
            case SearchlightMessenger::mtSLBalance:
            case SearchlightMessenger::mtSLMeta:
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
        auto qiter = served_queries_.emplace(query_id,
                std::make_shared<QueryContext>(query)).first;
        qiter->second->stats_.forwards_sent_.resize(query->getInstancesCount());

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

    // Stats
    std::ostringstream os;
    os << "Query statistics follows (id=" << query_id << "):\n";
    iter->second->stats_.print(os);
    logger->info(os.str(), LOG4CXX_LOCATION);

    // Unblock slots waiting for requests (with error)
    for (auto slot_id: iter->second->slots_used_) {
        auto &slot = slots_[slot_id];
        std::lock_guard<std::mutex> slot_lock{slot.mtx_};
        slot.error_ = true;
        slot.cond_.notify_one();
    }
    iter->second->slots_used_.clear();

    // Get rid of the context
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

void SearchlightMessenger::SetDistributedMapUpdateFrequency(
        const boost::shared_ptr<Query> &query,
        const std::string &array_name, int map_update_freq) {
    QueryContextPtr query_ctx = GetQueryContext(query->getQueryID(), true);
    auto array = GetRegisteredArray(query_ctx, array_name);
    array->map_broadcast_period_ = map_update_freq;
    LOG4CXX_DEBUG(logger, "Set update frequency: array=" << array_name
            << ", freq=" << map_update_freq);
}

SearchlightMessenger::QueryContext::ServedArrayPtr
SearchlightMessenger::GetRegisteredArray(
        const SearchlightMessenger::QueryContextPtr &query_ctx,
        const std::string &name) const {
    const auto &iter = query_ctx->reg_arrays_.find(name);
    if (iter == query_ctx->reg_arrays_.end()) {
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

void SearchlightMessenger::ReturnMessageSlot(size_t slot,
        QueryContextPtr &ctx) {
    std::lock_guard<std::mutex> lock(mtx_);
    assert(slots_[slot].used_);
    slots_[slot].used_ = false;
    slots_[slot].data_ready_ = false;
    slots_[slot].data_ = nullptr;
    slots_[slot].error_ = false;
    ctx->slots_used_.erase(slot);
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
    factory->addMessageType(mtSLControl,
            MessageCreator(&CreateMessageRecord),
            boost::bind(&SearchlightMessenger::HandleGeneralMessage, this, _1));
    factory->addMessageType(mtSLBalance,
            MessageCreator(&CreateMessageRecord),
            boost::bind(&SearchlightMessenger::HandleGeneralMessage, this, _1));
    factory->addMessageType(mtSLMeta,
            MessageCreator(&CreateMessageRecord),
            boost::bind(&SearchlightMessenger::HandleMetaMessage, this, _1));
}

void SearchlightMessenger::GetDistrChunksInfo(
        const boost::shared_ptr<Query> &query,
        const std::string &array_name,
        const CoordinateSet& chunks,
        std::vector<int> instance_counts) const {

    std::lock_guard<std::mutex> lock(mtx_);
    QueryContextPtr query_ctx = GetQueryContext(query->getQueryID(), false);
    try {
        QueryContext::ServedArrayPtr array =
                GetRegisteredArray(query_ctx, array_name);
        for (const auto &chunk: chunks) {
            const auto &inst_iter = array->chunks_map_.find(chunk);
            if (inst_iter != array->chunks_map_.end()) {
                for (auto inst: inst_iter->second) {
                    assert(inst < instance_counts.size());
                    ++instance_counts[inst];
                }
            }
        }
    } catch (const scidb::SystemException &ex) {
        // Array not found; just don't update counters
    }
}

bool SearchlightMessenger::RequestChunk(const boost::shared_ptr<Query> &query,
        InstanceID inst, const std::string &array_name, const Coordinates &pos,
        AttributeID attr, Chunk *chunk) {
    // check if the array is registered
    QueryContextPtr query_ctx = GetQueryContext(query->getQueryID(), true);
    GetRegisteredQuery(query_ctx); // just error checking
    QueryContext::ServedArrayPtr array =
            GetRegisteredArray(query_ctx, array_name);

    // message params
    const bool confirm_only = (chunk == nullptr);
    const size_t slot_num = GetMessageSlot();
    MessageSlot &msg_slot = slots_[slot_num];
    {
        std::lock_guard<std::mutex> lock{mtx_};
        query_ctx->slots_used_.insert(slot_num);
    }

    std::shared_ptr<ChunkRequestData> chunk_req =
            std::make_shared<ChunkRequestData>();
    chunk_req->chunk_ = chunk;
    chunk_req->array_ = array->array_;
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
    if (logger->isTraceEnabled()) {
        std::ostringstream os;
        os <<"Requesting chunk: qid=" << query->getQueryID() << ", name=" <<
                array_name << ", pos=" << pos << ", attr=" << attr <<
                ", data=" <<
                std::boolalpha << !confirm_only;
        logger->trace(os.str(), LOG4CXX_LOCATION);
    }

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(inst, msg);
    query_ctx->stats_.msgs_sent_++;

    // now we have to wait
    {
        const auto wait_start_time = std::chrono::steady_clock::now();

        std::unique_lock<std::mutex> slot_lock(msg_slot.mtx_);
        while (!msg_slot.data_ready_ && !msg_slot.error_) {
            msg_slot.cond_.wait(slot_lock);
        }

        const auto wait_end_time = std::chrono::steady_clock::now();
        query_ctx->stats_.total_wait_time_ +=
                std::chrono::duration_cast<std::chrono::microseconds>(
                        wait_end_time - wait_start_time);

        // Check for error
        if (msg_slot.error_) {
            ReturnMessageSlot(slot_num, query_ctx);
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR, SCIDB_LE_ILLEGAL_OPERATION)
                    << "Error during chunk request, possible query abort";
        }

        // If we requested data, the chunk is local now (empty or nor)
        if (!confirm_only) {
            array->chunks_map_[pos].push_back(query->getInstanceID());
            if (array->map_broadcast_period_) {
                array->last_fetched_chunks_.push_back(pos);
                if (array->last_fetched_chunks_.size() ==
                        array->map_broadcast_period_) {
                    BroadcastDistrUpdate(query, array_name,
                            array->last_fetched_chunks_);
                    array->last_fetched_chunks_.clear();
                }
            }
        }
        // safe to unlock since we write only once to the slot
    }

    // got the data
    const bool chunk_empty = chunk_req->empty_;
    ReturnMessageSlot(slot_num, query_ctx);

    return !chunk_empty;
}

void SearchlightMessenger::HandleSLChunkRequest(
        const boost::shared_ptr<MessageDescription> &msg_desc) {

    // some correctness checking
    const QueryContextPtr query_ctx =
            GetQueryContext(msg_desc->getQueryId(), true);
    const boost::shared_ptr<Query> query = GetRegisteredQuery(query_ctx);
    query_ctx->stats_.msgs_received_++;

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
    if (logger->isTraceEnabled()) {
        std::ostringstream os;
        os <<"Got chunk request: qid=" << query->getQueryID() << ", name=" <<
                array_name << ", pos=" << pos << ", attr=" << attr <<
                ", data=" <<
                std::boolalpha << !confirm_only;
        logger->trace(os.str(), LOG4CXX_LOCATION);
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
        if (array->getArrayDesc().getEmptyBitmapAttribute() &&
                !chunk.getAttributeDesc().isEmptyIndicator()) {
            empty_bitmap = chunk.getEmptyBitmap();
        }
        chunk.compress(*buffer, empty_bitmap);
        empty_bitmap.reset(); // apparently, reset() is mandatory here (bug?)

        reply->set_compression_method(buffer->getCompressionMethod());
        reply->set_decompressed_size(buffer->getDecompressedSize());
        reply->set_count(chunk.isCountKnown() ? chunk.count() : 0);

        // stats
        query_ctx->stats_.chunks_sent_++;
        query_ctx->stats_.chunk_data_sent_ += buffer->getSize();
    }

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->sendPhysical(msg_desc->getSourceInstanceID(), msg);
    query_ctx->stats_.msgs_sent_++;
}

void SearchlightMessenger::HandleSLChunk(
        const boost::shared_ptr<MessageDescription> &msg_desc) {

    // some correctness checking
    const QueryContextPtr query_ctx =
            GetQueryContext(msg_desc->getQueryId(), true);
    const boost::shared_ptr<Query> query = GetRegisteredQuery(query_ctx);
    query_ctx->stats_.msgs_received_++;

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
        query_ctx->stats_.chunks_received_++;
        query_ctx->stats_.chunk_data_received_ +=
                compressed_buffer->getSize();

        Chunk *chunk = chunk_req->chunk_;

        chunk->decompress(*compressed_buffer);
        chunk->setCount(record->count());
        /*
         *  Note, we do not call chunk->write() here, since the caller might
         *  want to do something else. write() might cause unpin() for some
         *  types of chunks (e.g., LRU Chunk) and a potential swap-out. This
         *  would be inefficient. So, write()/unpin() is the caller's
         *  responsibility.
         */
    }

    // log
    if (logger->isTraceEnabled()) {
        std::ostringstream os;
        os <<"Got chunk: qid=" << query->getQueryID() << ", name=" <<
                chunk_req->array_->getName() << ", empty=" <<
                std::boolalpha << chunk_req->empty_;
        logger->trace(os.str(), LOG4CXX_LOCATION);
    }

    std::lock_guard<std::mutex> lock(slot.mtx_);
    slot.data_ready_ = true;
    slot.cond_.notify_one();
}

void SearchlightMessenger::HandleMetaMessage(
        const boost::shared_ptr<MessageDescription> &msg_desc) {
    // some correctness checking
    const QueryContextPtr query_ctx =
            GetQueryContext(msg_desc->getQueryId(), true);
    query_ctx->stats_.msgs_received_++;

    // get the record
    boost::shared_ptr<SearchlightMeta> record =
            boost::dynamic_pointer_cast<SearchlightMeta>(msg_desc->getRecord());
    if (!record) {
        throw SYSTEM_EXCEPTION(scidb::SCIDB_SE_NETWORK,
                scidb::SCIDB_LE_INVALID_MESSAGE_FORMAT) <<
                        msg_desc->getMessageType();
    }
    const InstanceID source_inst = msg_desc->getSourceInstanceID();

    switch (record->type()) {
        case SearchlightMeta::DYNAMIC_DISTR: {
            std::lock_guard<std::mutex> lock{mtx_};
            const auto &distr_info = record->chunks();
            auto array = GetRegisteredArray(query_ctx, distr_info.array_name());
            for (const auto &chunk_pos: distr_info.chunks()) {
                Coordinates pos{chunk_pos.position().begin(),
                    chunk_pos.position().end()};
                array->chunks_map_[pos].push_back(source_inst);
            }
            break;
        }
    }
}

void SearchlightMessenger::HandleGeneralMessage(
        const boost::shared_ptr<MessageDescription> &msg_desc) {

    // get the query
    const QueryContextPtr query_ctx =
            GetQueryContext(msg_desc->getQueryId(), true);
    const boost::shared_ptr<Query> query = GetRegisteredQuery(query_ctx);
    query_ctx->stats_.msgs_received_++;

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
        LOG4CXX_TRACE(logger,
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
        const std::vector<int64_t> &vals,
        const std::vector<int64_t> &add_vals) const {

    // determine the coordinator
    const InstanceID coord_id = query->getCoordinatorID();
    if (coord_id == scidb::INVALID_INSTANCE) {
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

    // fill the record (values in mins, additional values in aux_info)
    VarAssignment &var_asgn = *record->mutable_solution();
    PackAssignment(vals, {}, var_asgn);
    if (!add_vals.empty()) {
        for (int64_t v: add_vals) {
            var_asgn.add_aux_info(v);
        }
    }

    // log
    LOG4CXX_TRACE(logger, "Sending solution to the coordinator: qid=" <<
            query->getQueryID());

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(coord_id, msg);

    // stats
    const QueryContextPtr query_ctx =
            GetQueryContext(query->getQueryID(), true);
    query_ctx->stats_.msgs_sent_++;
}

void SearchlightMessenger::ReportIdleSolver(
        const boost::shared_ptr<Query> &query, uint64_t solver_id) const {
    // determine the coordinator
    const InstanceID coord_id = query->getCoordinatorID();
    if (coord_id == scidb::INVALID_INSTANCE) {
        LOG4CXX_ERROR(logger, "Attempting to report idle coordinator");
        return;
    }

    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg(new SearchlightMessageDesc);
    msg->initRecord(mtSLControl);
    msg->setQueryID(query->getQueryID());
    boost::shared_ptr<SearchlightControl> record =
            msg->getRecord<SearchlightControl>();

    // fill the record
    record->set_type(SearchlightControl::SEARCH_IDLE);
    record->add_id(solver_id);

    // log
    LOG4CXX_DEBUG(logger, "Reporting idle search: qid=" << query->getQueryID());

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(coord_id, msg);

    // stats
    const QueryContextPtr query_ctx =
            GetQueryContext(query->getQueryID(), true);
    query_ctx->stats_.msgs_sent_++;
}

void SearchlightMessenger::ReportFinValidator(
        const boost::shared_ptr<Query> &query) const {
    // determine the coordinator
    const InstanceID coord_id = query->getCoordinatorID();
    if (coord_id == scidb::INVALID_INSTANCE) {
        LOG4CXX_ERROR(logger, "Attempting to report idle validator at "
                "the coordinator");
        return;
    }

    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg(new SearchlightMessageDesc);
    msg->initRecord(mtSLControl);
    msg->setQueryID(query->getQueryID());
    boost::shared_ptr<SearchlightControl> record =
            msg->getRecord<SearchlightControl>();

    // fill the record
    record->set_type(SearchlightControl::VALIDATOR_LOCAL_FIN);

    // log
    LOG4CXX_DEBUG(logger, "Reporting locally finished validator: qid=" <<
            query->getQueryID());

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(coord_id, msg);

    // stats
    const QueryContextPtr query_ctx =
            GetQueryContext(query->getQueryID(), true);
    query_ctx->stats_.msgs_sent_++;
}

void SearchlightMessenger::DispatchHelper(const boost::shared_ptr<Query> &query,
        uint64_t helper, uint64_t helpee, InstanceID dest) const {
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance> record =
            msg->getRecord<SearchlightBalance>();

    // fill the record
    record->set_type(SearchlightBalance::HELPER_DISPATCH);
    record->add_id(helpee);
    record->add_id(helper);

    // log
    LOG4CXX_DEBUG(logger, "Sending a helper: id=" <<
            helper << ", dest=" << dest);

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(dest, msg);
}

void SearchlightMessenger::RejectHelp(const boost::shared_ptr<Query> &query,
        const std::vector<InstanceID> &ids, uint64_t src, bool hard) const {
    // determine the coordinator
    const InstanceID coord_id = query->getCoordinatorID();
    if (coord_id == scidb::INVALID_INSTANCE) {
        LOG4CXX_ERROR(logger, "Attempting to reject help at the coordinator");
        return;
    }

    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance> record =
            msg->getRecord<SearchlightBalance>();

    // Fill the record
    record->set_type(hard ? SearchlightBalance::REJECT_HELPER_HARD :
            SearchlightBalance::REJECT_HELPER_SOFT);
    record->add_id(src);
    for (auto id: ids) {
        record->add_id(id);
    }

    // Log
    LOG4CXX_DEBUG(logger, "Rejecting help: hard=" << (hard ? "true" : "false"));

    // Send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(coord_id, msg);
}

void SearchlightMessenger::DispatchWork(const boost::shared_ptr<Query> &query,
        const CandidateVector &work,
        uint64_t solver, InstanceID dest) const {
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance> record =
            msg->getRecord<SearchlightBalance>();

    // fill the record
    record->set_type(SearchlightBalance::HELP_LOAD);
    FillBalanceMessage(*record, work);
    record->add_id(solver);

    // log
    LOG4CXX_DEBUG(logger, "Sending additional load to helper: qid=" <<
            query->getQueryID() << ", helper=" << solver);

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(dest, msg);
}

void SearchlightMessenger::RequestHelp(const boost::shared_ptr<Query> &query,
                 uint64_t solver) const {
    // prepare the message
    const InstanceID coord_id = query->getCoordinatorID();
    if (coord_id == scidb::INVALID_INSTANCE) {
        LOG4CXX_ERROR(logger, "Attempting to accept help at the coordinator");
        return;
    }
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance>  record =
            msg->getRecord<SearchlightBalance>();

    // Fill the record
    record->set_type(SearchlightBalance::REQ_HELP);
    record->add_id(solver);
    // log
    LOG4CXX_DEBUG(logger, "Sending request for help: " << solver);
    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(coord_id, msg);
}

void SearchlightMessenger::AcceptHelp(const boost::shared_ptr<Query> &query,
        uint64_t inst) const {
    // prepare the message
    const InstanceID coord_id = query->getCoordinatorID();
    if (coord_id == scidb::INVALID_INSTANCE) {
        LOG4CXX_ERROR(logger, "Attempting to accept help at the coordinator");
        return;
    }

    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance>  record =
            msg->getRecord<SearchlightBalance>();

    // Fill the record
    record->set_type(SearchlightBalance::ACCEPT_HELP);
    record->add_id(inst);

    // log
    LOG4CXX_DEBUG(logger, "Sending help accept: helper=" << inst);

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(coord_id, msg);
}

void SearchlightMessenger::ForwardCandidates(
        const boost::shared_ptr<Query> &query,
        const CandidateVector &cands,
        const std::vector<int64_t> &zones,
        InstanceID dest,
        uint64_t forw_id) const {
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance> record =
            msg->getRecord<SearchlightBalance>();

    // fill the record
    record->set_type(SearchlightBalance::CANDIDATE_FORWARD);
    FillBalanceMessage(*record, cands, zones);
    record->add_id(forw_id);
    {
        std::lock_guard<std::mutex> lock{mtx_};
        GetQueryContext(query->getQueryID(), false)->
                stats_.forwards_sent_[dest] += cands.size();
    }

    // log
    LOG4CXX_TRACE(logger, "Forwarding candidates: dest=" << dest
            << ", forw_id=" << forw_id);

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(dest, msg);
}

void SearchlightMessenger::BroadcastDistrUpdate(
        const boost::shared_ptr<Query> &query,
        const std::string &array_name,
        const std::vector<Coordinates> &chunks) const {
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLMeta);
    boost::shared_ptr<SearchlightMeta> record =
            msg->getRecord<SearchlightMeta>();

    // fill it
    record->set_type(SearchlightMeta::DYNAMIC_DISTR);
    SearchlightMeta::ChunksInfo *map = record->mutable_chunks();
    map->set_array_name(array_name);
    for (const auto &chunk_pos: chunks) {
        auto chunk = map->add_chunks();
        for (auto pos: chunk_pos) {
            chunk->add_position(pos);
        }
    }

    // log
    LOG4CXX_DEBUG(logger, "Broadcasting distribution update...");

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->broadcastLogical(msg);
}

void SearchlightMessenger::SendBalanceResult(
        const boost::shared_ptr<Query> &query,
        InstanceID dest, uint64_t forw_id, bool result,
		const std::vector<int64_t> &add_vals) const {

    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance> record =
            msg->getRecord<SearchlightBalance>();

    // fill the record
    record->set_type(SearchlightBalance::BALANCE_RESULT);
    record->add_id(forw_id);
    record->set_result(result);
    if (!add_vals.empty()) {
		VarAssignment *var_asgn = record->add_load();
		PackAssignment(add_vals, {}, *var_asgn);
    }

    // log
    LOG4CXX_TRACE(logger, "Sending balancing result: dest=" << dest
            << ", id=" << forw_id << ", result=" << result);

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->send(dest, msg);
}

void SearchlightMessenger::BroadcastFinishSearch(
        const boost::shared_ptr<Query> &query) const {
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg(new SearchlightMessageDesc);
    msg->initRecord(mtSLControl);
    msg->setQueryID(query->getQueryID());
    boost::shared_ptr<SearchlightControl> record =
            msg->getRecord<SearchlightControl>();

    // fill the record
    record->set_type(SearchlightControl::END_SEARCH);

    // log
    LOG4CXX_DEBUG(logger, "Broadcasting end-of-search: qid=" <<
            query->getQueryID());

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->broadcastLogical(msg);

    // stats
    const QueryContextPtr query_ctx =
            GetQueryContext(query->getQueryID(), true);
    query_ctx->stats_.msgs_sent_++;
}

void SearchlightMessenger::BroadcastValidatorInfo(
        const boost::shared_ptr<Query> &query,
        size_t cands_num) const {
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance> record =
            msg->getRecord<SearchlightBalance>();

    // fill the record
    record->set_type(SearchlightBalance::VALIDATOR_INFO);
    record->add_id(cands_num);

    // log
    LOG4CXX_TRACE(logger, "Broadcasting validator info: cands=" << cands_num);

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->broadcastLogical(msg);
}

void SearchlightMessenger::BroadcastRD(const boost::shared_ptr<Query> &query,
                                       double rd) const {
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg =
            PrepareMessage(query->getQueryID(), mtSLBalance);
    boost::shared_ptr<SearchlightBalance> record =
            msg->getRecord<SearchlightBalance>();

    // fill the record
    record->set_type(SearchlightBalance::LRD);
    record->set_lrd(rd);

    // log
    LOG4CXX_TRACE(logger, "Broadcasting new RD: rd=" << rd);

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->broadcastLogical(msg);
}

void SearchlightMessenger::BroadcastCommit(
        const boost::shared_ptr<Query> &query) const {
    // prepare the message
    boost::shared_ptr<scidb::MessageDesc> msg(new SearchlightMessageDesc);
    msg->initRecord(mtSLControl);
    msg->setQueryID(query->getQueryID());
    boost::shared_ptr<SearchlightControl> record =
            msg->getRecord<SearchlightControl>();

    // fill the record
    record->set_type(SearchlightControl::COMMIT);

    // log
    LOG4CXX_DEBUG(logger, "Broadcasting search commit: qid=" <<
            query->getQueryID());

    // send
    NetworkManager *network_manager = NetworkManager::getInstance();
    network_manager->broadcastLogical(msg);

    // stats
    const QueryContextPtr query_ctx =
            GetQueryContext(query->getQueryID(), true);
    query_ctx->stats_.msgs_sent_++;
}

void SearchlightMessenger::Synchronize(const boost::shared_ptr<Query> &query) {
    LOG4CXX_DEBUG(logger, "Requesting SL Messenger synchronization...");
    // The same id should be okay here, even in case of subsequent barriers
    scidb::syncBarrier(0, query);
    LOG4CXX_DEBUG(logger, "Completed SL Messenger synchronization...");
}

void SearchlightMessenger::PackAssignment(
        const std::vector<int64_t> &mins, const std::vector<int64_t> &maxs,
        VarAssignment &msg) {
    const bool is_range = !maxs.empty();
    assert(!is_range || mins.size() == maxs.size());
    for (size_t i = 0; i < mins.size(); i++) {
        msg.add_var_min(mins[i]);
        if (is_range) {
            msg.add_var_max(maxs[i]);
        }
    }
}

void SearchlightMessenger::UnpackAssignment(const VarAssignment &msg,
        std::vector<int64_t> &mins, std::vector<int64_t> &maxs) {
    const bool is_range = msg.var_max_size();
    mins.resize(msg.var_min_size());
    if (is_range) {
        assert(msg.var_min_size() == msg.var_max_size());
        maxs.resize(msg.var_max_size());
    }
    for (int i = 0; i < msg.var_min_size(); i++) {
        mins[i] = msg.var_min(i);
        if (is_range) {
            maxs[i] = msg.var_max(i);
        }
    }
}

void SearchlightMessenger::UnpackAssignment(const VarAssignment &msg,
        CandidateAssignment &asgn) {
    UnpackAssignment(msg, asgn.var_asgn_.mins_, asgn.var_asgn_.maxs_);
    if (msg.rel_const_size()) {
        asgn.relaxed_constrs_.resize(msg.rel_const_size());
        for (int i = 0; i < msg.rel_const_size(); ++i) {
            asgn.relaxed_constrs_[i] = msg.rel_const(i);
        }
        asgn.best_rd_ = msg.rd();
    }
}

void SearchlightMessenger::PackAssignment(const CandidateAssignment &asgn,
        VarAssignment &msg, const std::vector<int64_t> &aux) {
    const LiteVarAssignment &var_asgn = asgn.var_asgn_;
    PackAssignment(var_asgn.mins_, var_asgn.maxs_, msg);
    // Aux info, if any
    for (const auto x: aux) {
        msg.add_aux_info(x);
    }
    // Relaxed constraint info
    for (const auto x: asgn.relaxed_constrs_) {
        msg.add_rel_const(x);
    }
    // Best RD
    msg.set_rd(asgn.best_rd_);
}

void SearchlightMessenger::FillBalanceMessage(SearchlightBalance &msg,
        const CandidateVector &asgns,
        const std::vector<int64_t> &aux) {
    if (aux.empty()) {
        FillBalanceMessage(msg, asgns);
    } else {
		std::vector<int64_t> aux_info(1);
		for (size_t i = 0; i < asgns.size(); ++i) {
			VarAssignment *var_asgn = msg.add_load();
			aux_info[0] = aux[i];
			PackAssignment(asgns[i], *var_asgn, aux_info);
		}
    }
}

void SearchlightMessenger::FillBalanceMessage(SearchlightBalance &msg,
        const CandidateVector &asgns) {
    for (size_t i = 0; i < asgns.size(); ++i) {
        VarAssignment *var_asgn = msg.add_load();
        PackAssignment(asgns[i], *var_asgn, {});
    }
}

boost::shared_ptr<scidb::MessageDesc> SearchlightMessenger::PrepareMessage(
        QueryID qid,
        SearchlightMessageType type) const {
    // Create message
    boost::shared_ptr<scidb::MessageDesc> msg(new SearchlightMessageDesc);
    msg->initRecord(type);
    msg->setQueryID(qid);

    // Update stats
    const QueryContextPtr query_ctx = GetQueryContext(qid, true);
    query_ctx->stats_.msgs_sent_++;

    return msg;
}
}
