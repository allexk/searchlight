/*
 * Copyright 2016, Brown University, Providence, RI.
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
 * @file synopsis.cpp
 * The implementation of the basic synopsis stuff
 *
 * @author Alexander Kalinin
 */

#include "synopsis.h"

namespace searchlight {
void AggCellItemReader::InitIterators() {
    /*
     * Here we assume the attributes have been discovered and checked.
     */
    assert(attributes_.count("min") && attributes_.count("max") &&
           attributes_.count("sum") && attributes_.count("count"));
    /*
     * One thing to consider here. Creating item iterator results in
     * fetching the first array chunk, which might create a small performance
     * penalty. Another solution is to use array iterators and create
     * chunk iterators when fetching a synopsis cell.
     */
    min_it_   = array_->getItemIterator(attributes_.find("min")->second);
    max_it_   = array_->getItemIterator(attributes_.find("max")->second);
    sum_it_   = array_->getItemIterator(attributes_.find("sum")->second);
    count_it_ = array_->getItemIterator(attributes_.find("count")->second);
}

void AggCellItemReader::Read(const Coordinates &pos, AggCell &cell) {
    // Init iterators (first time only)
    if (!count_it_) {
        InitIterators();
    }
    if (!count_it_->setPosition(pos) || count_it_->isEmpty()) {
        // No chunk in the array -- assume the cell is empty
        cell.count_ = 0;
    } else {
        // should be a non-empty chunk for sure
        if (!min_it_->setPosition(pos) || !max_it_->setPosition(pos) ||
                !sum_it_->setPosition(pos)) {
            std::ostringstream err_msg;
            err_msg << "Cannot get info from synopsis, pos=(";
            std::copy(pos.begin(), pos.end(),
                    std::ostream_iterator<Coordinate>(err_msg, ", "));
            err_msg << ")";

            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
        }
        // Read the cell
        cell.count_ = count_it_->getItem().getUint64();
        /*
         * We might still have an empty cell, e.g., when all values are null
         * there.
         */
        if (cell.count_) {
            cell.min_ = min_it_->getItem().getDouble();
            cell.max_ = max_it_->getItem().getDouble();
            cell.sum_ = sum_it_->getItem().getDouble();
            if (std::isnan(cell.min_) || std::isnan(cell.max_) ||
                    std::isnan(cell.sum_)) {
                // Treat cells with NANs as "corrupted"; just count as empty
                cell.count_ = 0;
            }
        }
    }
}

void SeqCellItemReader::InitIterators() {
    /*
     * Here we assume the attributes have been discovered and checked.
     */
    assert(attributes_.count("low") && attributes_.count("high"));
    /*
     * One thing to consider here. Creating item iterator results in
     * fetching the first array chunk, which might create a small performance
     * penalty. Another solution is to use array iterators and create
     * chunk iterators when fetching a synopsis cell.
     */
    low_it_  = array_->getItemIterator(attributes_.find("low")->second);
    high_it_ = array_->getItemIterator(attributes_.find("high")->second);
    /*
     * The number of DFT components depends on the last dimension size, which
     * is used to store them.
     */
    const auto &dims = array_->getArrayDesc().getDimensions();
    coords_num_ = dims.back().getLength();
    read_pos_.resize(dims.size());
}

void SeqCellItemReader::Read(const Coordinates &pos, SeqCell &cell) {
    // Init iterators (first time only)
    if (!low_it_) {
        InitIterators();
    }
    // A small check that the number of dimensions is right
    assert(array_->getArrayDesc().getDimensions().size() == pos.size() + 1 &&
           read_pos_.size() == pos.size() + 1);
    std::memcpy(read_pos_.data(), pos.data(), sizeof(Coordinate) * pos.size());
    for (size_t i = 0; i < coords_num_; ++i) {
        read_pos_.back() = i;
        // All coordinates should be present
        if (low_it_->setPosition(read_pos_) && !low_it_->isEmpty() &&
                high_it_->setPosition(read_pos_) && !high_it_->isEmpty()) {
            if (i == 0) {
                cell.mbr_.low_.resize(coords_num_);
                cell.mbr_.high_.resize(coords_num_);
            }
            cell.mbr_.low_[i]  = low_it_->getItem().getDouble();
            cell.mbr_.high_[i] = high_it_->getItem().getDouble();
        } else if (i == 0) {
            // No components: consider the cell empty
            break;
        } else {
            // One of the components is absent? Indicate an error.
            assert(i > 0);
            std::ostringstream err_msg;
            err_msg << "Cannot get info from synopsis, pos=(";
            std::copy(read_pos_.begin(), read_pos_.end(),
                    std::ostream_iterator<Coordinate>(err_msg, ", "));
            err_msg << ")";
            throw SYSTEM_EXCEPTION(SCIDB_SE_OPERATOR,
                    SCIDB_LE_ILLEGAL_OPERATION) << err_msg.str();
        }

    }
}
} /* namespace searchlight */
