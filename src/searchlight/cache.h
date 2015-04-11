/*
 * Copyright 2013, Brown University, Providence, RI.
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
 * @file cache.h
 * This file contains the definition of a cache structure that can be used to
 * cache objects between queries.
 *
 * The cache is thread-safe. The typical usage is to instantiate it as a
 * static in each class that wants to cache something.
 *
 * The cache always contains shared pointers instead of in-place values.
 * This is done for the purpose of sharing values between clients as well
 * as safe clearing. When Clear() is called, values for which external
 * pointers exist remain live until they are still in use.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_CACHE_H_
#define SEARCHLIGHT_CACHE_H_

#include <mutex>
#include <unordered_map>

namespace searchlight {

/**
 * Cache, storing objects of type Val with keys of type Key.
 *
 * Thread-safe. Val must be copy-constructible
 * and have a default constructor (for default value). The best way to use the
 * class is store pointers or (better) RAII pointers, like std::shared_ptr.
 */
template <typename Key, typename Val>
class SearchlightCache {
private:
    typedef std::unordered_map<Key, std::shared_ptr<Val>> CacheType;
    std::mutex mtx_;
    CacheType cache_;

public:
    /**
     * Iterator type for the cache.
     */
    typedef typename CacheType::iterator iterator;

    /**
     * Value with the shared pointer.
     */
    typedef typename std::shared_ptr<Val> SharedValue;

    /**
     * Check if the object is in the cache.
     *
     * @param key search key for the object
     * @return true, if the object exists; false, otherwise
     */
    bool CheckCache(const Key &key) const {
        std::lock_guard<std::mutex> lock{mtx_};
        return cache_.find(key) != cache_.end();
    }

    /**
     * Puts a new object in the map.
     *
     * It either puts a new object or returns the existing element. The
     * result is reflected in the return value.
     *
     * If the object already exists, args remain untouched and no
     * object is constructed.
     *
     * Note, this function returns the shared value to avoid problems with
     * concurrent updates. If a concurrent update is made, the iterator or
     * a reference to the value might get invalidated.
     *
     * @param args key and any remaining value constructor arguments
     * @return shared value and indicator if the insert took place
     */
    template <typename... Args>
    std::pair<SharedValue, bool> emplace(const Key &key, Args&&... args) {
        std::lock_guard<std::mutex> lock{mtx_};
        auto iter = cache_.find(key);
        if (iter == cache_.end()) {
            auto emp_iter = cache_.emplace(std::piecewise_construct,
                    std::forward_as_tuple(key),
                    std::forward_as_tuple(std::make_shared<Val>(
                            std::forward<Args>(args)...)));
            return std::make_pair(emp_iter.first->second, emp_iter.second);
        } else {
            return std::make_pair(iter->second, false);
        }
    }

    /**
     * Remove object with the given key.
     *
     * If the object is not in the cache, do nothing.
     *
     * @param key object search key
     */
    void Remove(const Key &key) {
        std::lock_guard<std::mutex> lock{mtx_};
        cache_.erase(key);
    }

    /**
     * Clear the cache.
     */
    void Clear() {
        std::lock_guard<std::mutex> lock{mtx_};
        cache_.clear();
    }
};

} /* namespace searchlight */
#endif /* SEARCHLIGHT_CACHE_H_ */
