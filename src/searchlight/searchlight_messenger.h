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
 * @file searchlight_messenger.h
 *
 * This file describes the Searchlight messenger. This class is responsible
 * with communicating with other instances to request chunks, send/receive
 * control messages, etc.
 *
 * The class uses internal SciDb NetworkManages services, including
 * networking and message formats. It registers new UDTs for messages, so
 * extra-care should be taken to avoid clashes with other message types,
 * e.g., MPI ones.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SEARCHLIGHT_MESSENGER_H_
#define SEARCHLIGHT_SEARCHLIGHT_MESSENGER_H_

#include "base.h"
#include "scidb_inc.h"

#include <mutex>
#include <condition_variable>
#include <unordered_map>

namespace searchlight {

class VarAssignment;
class SearchlightBalance;

/**
 * This class represents SearchlightMessenger. It allows communticating with
 * other instances running Searchlight to exchange chunks, control
 * messages, etc. This class is a singleton instantiated by the first user.
 * Then, each user should call RegisterQuery() to initialize it (only the first
 * call for the given query will come through).
 * The Messenger will deactivate itself automatically via the query finalizer.
 *
 * Messenger allows exchanging chunks of the arrays. For this, each array
 * must be registered at the Messenger. Otherwise, a request will cause an
 * exception.
 *
 * We will reuse Fetch and Chunk protobuf messages, already implemented at
 * SciDb. Unfortunately, we need to define our own handler and creator
 * callbacks, since mtFetch/mtChunk messages are considered system and, thus,
 * handled by SciDb internally. We re-use protobuf messages as typed data
 * carriers only.
 */
class SearchlightMessenger : public scidb::Singleton<SearchlightMessenger> {
public:

    /**
     * Searchlight message types.
     *
     * We start them with 50 to avoid clashes with SciDb message types,
     * especially for future versions.
     */
    enum SearchlightMessageType {
        /** Chunk request from another instance */
        mtSLChunkRequest = 50,
        /** Chunk arrived from another instance */
        mtSLChunk,
        /** A solution arrived from another instance */
        mtSLSolution,
        /** Control message between instances and coordinator */
        mtSLControl,
        /** Balance message between instances */
        mtSLBalance
    };

    /**
     * Handler for a message. Handlers are called for message types they
     * are registered for. They make sense only for asynchronous messages,
     * like SLResult.
     *
     * @param first_argument source instance id
     * @param second_argument pointer to the message record
     */
    typedef std::function<
            void(InstanceID, const google::protobuf::Message *)>
        UserMessageHandler;

    /**
     *  Constructs the messenger. For now it just registers message handlers.
     */
    SearchlightMessenger() {
        RegisterHandlers(true);
        /*
         * We do not remove message handlers on destruction. First of all,
         * SciDb does not have API for that -- handlers remain until
         * termination. Secondly, we just throw an exception if a request came
         * out of context (query). That will happen during query validation.
         *
         * An alternative would be to establish empty handlers, so SciDb
         * would drop messages by itself. However, we won't detect errors
         * that way.
         */
    }

    /**
     * Initializes the messenger for a new query. Can be called several times
     * within the same query -- only the first call counts. The messenger
     * will deactivate itself automatically, when the query is finished.
     *
     * It is assumed that all arrays are registered before queries begin
     * flying around. Since requests check if an array is registered,
     * violating this assumption might cause data races. Moreover, such
     * behavior might cause unexpected bugs anyway, if no synchronization
     * between the instances is performed after the registration.
     *
     * @param query current query
     */
    void RegisterQuery(const boost::shared_ptr<Query> &query);

    /**
     * Registers array for chunk exchange. Instances can exchange chunks only
     * for arrays registered at the messenger. All arrays are unregistered
     * when the current query finishes.
     *
     * @param query caller's query
     * @param array array to register
     */
    void RegisterArray(const boost::shared_ptr<Query> &query,
            const ArrayPtr &array);

    /**
     * Requests a chunk from the the specified array. The caller has to
     * explicitly specify what instance to request the chunk from. If the
     * caller does not specify a pointer for the chunk to be written, only
     * the existence of the chunk at the specified position is confirmed.
     * If the chunk is specified, it is supposed to be pinned in memory.
     * This function does not un-pin it; that is the caller's responsibility.
     *
     * @param query caller's query
     * @param inst  instance to requests the chunk from (logical ID)
     * @param array_name name of the array
     * @param pos   position of the hunk
     * @param attr  chunk attribute
     * @param chunk pointer to the destination (nullptr means confirm-only)
     * @return true, if the chunk exists; false, no chunk found
     */
    bool RequestChunk(const boost::shared_ptr<Query> &query,
            InstanceID inst, const std::string &array_name,
            const Coordinates &pos, AttributeID attr, Chunk *chunk);

    /**
     * Reports the main solver (search process) as idle.
     *
     * Idle means the solver finished its portion of the search space and can
     * accept additional load.
     *
     * @param query the current query
     */
    void ReportIdleSolver(const boost::shared_ptr<Query> &query) const;

    /**
     * Reports the local validator that already finished its local job.
     *
     * @param query the current query
     */
    void ReportFinValidator(const boost::shared_ptr<Query> &query) const;

    /**
     * Synchronizes all Searchlight Messenger instances. This means it
     * creates a barrier across instances participating in the query.
     *
     * Current implementation uses Scatter/Gather's barrier to achieve this,
     * since it is exactly what we want, and sg() should not run concurrently
     * with any operators using the messenger.
     *
     * Caveat: only one thread per instance should call this function.
     * Otherwise, the result is undefined. One good point to do this is during
     * physical plan preparation or execute(), which SciDb runs in a single
     * thread.
     *
     * @param query the query context
     */
    static void Synchronize(const boost::shared_ptr<Query> &query);

    /**
     * Register handler for the specified message type. The handler will be
     * called when a message of this type arrives at this instance.
     *
     * @param query the current query
     * @param msg_type type of message to register the handler for
     * @param handler the handler to call on arrival
     */
    void RegisterUserMessageHandler(
            const boost::shared_ptr<Query> &query,
            SearchlightMessageType msg_type,
            const UserMessageHandler &handler);

    /**
     * Sends end-of-search control message to all instances.
     *
     * @param query the current query
     */
    void BroadcastFinishSearch(const boost::shared_ptr<Query> &query) const;

    /**
     * Sends search commit control message to all instances.
     *
     * @param query the current query
     */
    void BroadcastCommit(const boost::shared_ptr<Query> &query) const;

    /**
     * Sends a Searchlight solution to the coordinator. If eor (end-of-result)
     * is true, var_mins and var_maxs are ignored. There should be an
     * appropriate handler registered at the coordinator's Messenger instance.
     *
     * @param query current query
     * @param eor true, if the result has ended; false, otherwise
     * @param vals solution values
     */
    void SendSolution(const boost::shared_ptr<Query> &query,
            const std::vector<int64_t> &vals) const;

    /**
     * Dispatches a helper to a remote instance.
     *
     * @param query caller's query context
     * @param helper instance willing to help
     * @param dest instance in need of help
     */
    void DispatchHelper(const boost::shared_ptr<Query> &query,
            InstanceID helper, InstanceID dest) const;

    /**
     * Dispatches work for a remote helper.
     *
     * @param caller's query
     * @param work assignments to send
     * @param solver id of the helper solver
     */
    void DispatchWork(const boost::shared_ptr<Query> &query,
            const LiteAssignmentVector &work,
            InstanceID solver) const;

    /**
     * Sends a message rejecting help to the coordinator.
     *
     * @param query caller's query context
     * @param ids ids of helpers
     * @param hard true, if it's a hard reject; false, if soft
     */
    void RejectHelp(const boost::shared_ptr<Query> &query,
            const std::vector<InstanceID> &ids, bool hard) const;

    /**
     * Forwards a candidate solution to another validator.
     *
     * The idea behind the id of the forward is that after the validation
     * the validator will respond with the same id.
     *
     * @param query caller's query context
     * @param cands candidates to forward
     * @param dest destination validator
     * @param forw_id id of this forward
     */
    void ForwardCandidates(const boost::shared_ptr<Query> &query,
            const LiteAssignmentVector &cands,
            InstanceID dest,
            int forw_id) const;

    /**
     * Sends result of the balancing.
     *
     * @param query caller's query id
     * @param dest destination Searchlight
     * @param id id of the balancing load
     * @param result result status of the load
     */
    void SendBalanceResult(const boost::shared_ptr<Query> &query,
            InstanceID dest, int id, bool result) const;

    /**
     * Checks if a chunk has been fetched and stored locally.
     *
     * This function only checks for fetched remote chunks. It does not
     * check for chunks that are local due to their distribution
     * (i.e., the instance owns it).
     *
     * If the array has not been registered with this messenger, the result
     * is always false, which is consistent with the knowledge this messenger
     * possess.
     *
     * @param query caller's query context
     * @param array_name the name of the array
     * @param pos chunk position
     * @return true, if the chunk has been fetched; false, otherwise
     */
    bool ChunkFetched(const boost::shared_ptr<Query> &query,
            const std::string &array_name, const Coordinates &pos) const;

    /**
     * Packs min/max assignment into a message.
     *
     * The suitable message must already be provided. maxs might be nullptr
     * for point assignments (e.g., where min[i] == max[i]).
     *
     * @param mins minimum values
     * @param maxs maximum values (nullptr if not required)
     * @param msg message to store the values into
     */
    static void PackAssignment(const std::vector<int64_t> &mins,
            const std::vector<int64_t> *maxs, VarAssignment &msg);

    /**
     * Packs an assignment into a message.
     *
     * If the assignment's maxs is empty, only min values are taken from it.
     * This is useful for complete assignments (i.e., without ranges).
     *
     * @param asgn assignment to pack
     * @param msg message to store the values into
     */
    static void PackAssignment(const LiteVarAssignment &asgn,
            VarAssignment &msg);

    /**
     * Unpacks a message assignment into min/max vectors.
     *
     * If the caller does not provide maxs, then max values are ignored. This
     * is useful for point assignments.
     *
     * @param msg message to get the values from
     * @param mins vector for minimum values
     * @param maxs vector for maximum values (nullptr if not required)
     */
    static void UnpackAssignment(const VarAssignment &msg,
            std::vector<int64_t> &mins, std::vector<int64_t> *maxs);

    /**
     * Unpacks a message assignment into assignment.
     *
     * If the message does not contain maxs, only mins are unpacked. This is
     * useful for complete (non-range) assignments.
     *
     * @param msg message to get the values from
     * @param asgn Assignment to unpack to
     */
    static void UnpackAssignment(const VarAssignment &msg,
            LiteVarAssignment &asgn);

private:
    /*
     * Information about a chunk request. Store in a mesage slot.
     */
    struct ChunkRequestData {
        ArrayPtr array_;  // array
        Chunk *chunk_;    // chunk data (nullptr -- no data required)
        bool empty_;      // is the chunk empty?
    };

    /*
     * Message slot. Used by a requester to wait for the result and to return
     * data. Each caller gets a slot.
     */
    struct MessageSlot {
        bool used_ = true; // is slot being used? for reusing the slot

        // sync waiting
        std::mutex mtx_;
        std::condition_variable cond_;

        // response data (type depends on the caller, e.g., ChunkRequestData)
        void *data_ = nullptr;
        bool data_ready_ = false; // is data ready?
    };
    std::deque<MessageSlot> slots_;

    /*
     * Every query served by the messenger has its own context. A query must
     * be explicitly registered with the messenger via the register() function.
     * Subsequent RegisterQuery() calls for the same query will be ignored.
     * When the query finished its context is de-allocated via a finalizer.
     */
    struct QueryContext {
        // Represents an array with cached iterators
        struct ServedArray {
            // The array itself
            const ArrayPtr array_;

            // Cached iterators
            std::unordered_map<AttributeID,
                boost::shared_ptr<ConstArrayIterator>> iters_;

            // Chunks retrieved from other nodes
            std::unordered_set<Coordinates, CoordinatesHash> fetched_chunks_;

            ServedArray(const ArrayPtr &array) : array_(array) {}
        };
        typedef std::shared_ptr<ServedArray> ServedArrayPtr;

        // Statistics about all querie's requests
        struct Statistics {
            std::atomic<uint32_t> msgs_sent_{0};
            std::atomic<uint32_t> msgs_received_{0};

            std::atomic<uint32_t> chunks_sent_{0};
            std::atomic<uint32_t> chunks_received_{0};

            std::atomic<uint64_t> chunk_data_sent_{0};
            std::atomic<uint64_t> chunk_data_received_{0};

            std::chrono::microseconds total_wait_time_;

            void print(std::ostream &os) const {
                os << "\n\tMessages sent=" << msgs_sent_;
                os << "\n\tMessages received=" << msgs_received_;
                os << "\n\tChunks sent=" << chunks_sent_;
                os << "\n\tChunks received=" << chunks_received_;
                os << "\n\tChunk data sent=" << chunk_data_sent_;
                os << "\n\tChunk data received=" << chunk_data_received_;
                os << '\n';
            }
        };
        Statistics stats_; // Statistics for the query

        // Reference to the query
        const boost::weak_ptr<Query> query_;

        // Arrays registered for chunk exchange
        std::unordered_map<std::string, ServedArrayPtr> reg_arrays_;

        // user handlers for some messages
        std::unordered_map<MessageID, UserMessageHandler>
            message_handlers_;

        QueryContext(const boost::shared_ptr<Query> &query) :
            query_(query) {}
    };
    typedef std::shared_ptr<QueryContext> QueryContextPtr;
    std::unordered_map<QueryID, QueryContextPtr> served_queries_;

    // To guard the structure
    mutable std::mutex mtx_;

    // Chunk request handler
    void HandleSLChunkRequest(
            const boost::shared_ptr<MessageDescription> &msg_desc);

    // Chunk response handler
    void HandleSLChunk(
            const boost::shared_ptr<MessageDescription> &msg_desc);

    // Handler for a message that should be handled by a user handler
    void HandleGeneralMessage(
            const boost::shared_ptr<MessageDescription> &msg_desc);

    // Helper to fill in a balance message
    static void FillBalanceMessage(SearchlightBalance &msg,
            const LiteAssignmentVector &asgns);

    // Deactivates the messenger from the query
    void deactivate(const boost::shared_ptr<Query> &query);

    // Register message handlers. init = true means we are initializing.
    void RegisterHandlers(bool init);

    // Returns a new message slot (or re-uses an old one)
    uint32_t GetMessageSlot();

    // Returns a message slot
    void ReturnMessageSlot(size_t slot);

    // Helper to create a message and update stats
    boost::shared_ptr<scidb::MessageDesc> PrepareMessage(
            QueryID qid,
            SearchlightMessageType type) const;

    // Returns pointer to the registered array with name
    QueryContext::ServedArrayPtr GetRegisteredArray(
            const QueryContextPtr &query_ctx,
            const std::string &name) const;

    // Returns valid pointer to the registered query
    boost::shared_ptr<Query> GetRegisteredQuery(
            const QueryContextPtr &query_ctx) const;

    // Returns query context for the given id; if lock == true, grab the mutex
    QueryContextPtr GetQueryContext(QueryID query_id, bool lock) const;
};
}
#endif /* SEARCHLIGHT_SEARCHLIGHT_MESSENGER_H_ */
