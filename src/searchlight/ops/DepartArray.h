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
 * @file DepartArray.h
 *
 * This file describes "departitioning" array that is returned from the
 * depart() operator. The goal of this array is to ensure global view of the
 * underlying array. It does so by either fetching chunks from the local
 * array partition or requesting chunks from other instances. It uses
 * Searchlight Messenger for communication.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_DEPART_ARRAY_H_
#define SEARCHLIGHT_DEPART_ARRAY_H_

#include "../scidb_inc.h"
#include "../searchlight_messenger.h"

#include <unordered_set>

namespace searchlight {

/**
 * DepartArray that ensures global view of the underlying input. The schema
 * remains the same. The array fetches chunks situated on other instances and
 * brings them to the local instance. Chunks are cached by MemArray that
 * automatically supports LRU.
 *
 */
class DepartArray : public DelegateArray {
public:
    /**
     * Creates a new departitioning array.
     *
     * @param desc          this array descriptor
     * @param input         input array
     * @param input_distr   input array's distribution
     * @param query         current query
     */
    DepartArray(const ArrayDesc &desc, const ArrayPtr &input,
            const ArrayDistribution &input_distr,
            const boost::shared_ptr<Query> &query) :
        DelegateArray(desc, input, false),
        messenger_(SearchlightMessenger::getInstance()),
        cache_(desc, query),
        query_(query),
        input_distr_(input_distr) {

        // register array at the messenger
        messenger_->RegisterArray(query, input);
    }

    /**
     * Creates a new const departitioning array iterator.
     *
     * @param attr attribute id
     * @return const iterator for the specified attribute
     */
    virtual boost::shared_ptr<ConstArrayIterator>
        getConstIterator(AttributeID attr) const override;

private:
    // Cache part of the array
    struct CacheMemArray {
        // Cache for remote chunks
        MemArray array_;

        // Remote chunks we've seen (true means the chunk is non-empty)
        std::unordered_map<Address, bool, AddressHash> remote_chunks_;

        // Current requests
        std::unordered_set<Address, AddressHash> current_requests_;

        // For synchronization
        std::mutex mtx_;
        std::condition_variable cond_;

        // Init the cache via the provided schema and context
        CacheMemArray(const ArrayDesc &array,
                const boost::shared_ptr<Query> &query) :
                    array_(array, query) {}
    };

    // To access cache and other info
    friend class DepartArrayIterator;

    // Messenger for communication
    SearchlightMessenger * const messenger_;

    // Cache for fetched chunks
    mutable CacheMemArray cache_;

    // The current query
    const boost::weak_ptr<Query> query_;

    // Input distribution
    const ArrayDistribution input_distr_;
};

}

#endif /* SEARCHLIGHT_DEPART_ARRAY_H_ */
