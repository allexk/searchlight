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
 * @file scidb_inc.h
 * This is common file for other sources using SciDb. It contains
 * common definitions and includes from SciDb.
 *
 * @author Alexander Kalinin
 */

#ifndef SEARCHLIGHT_SCIDB_INC_H_
#define SEARCHLIGHT_SCIDB_INC_H_

#include <array/Metadata.h>
#include <array/Array.h>

using scidb::SystemCatalog;
using scidb::Config;
using scidb::ArrayDesc;
using scidb::Array;
using scidb::Attributes;
using scidb::AttributeID;
using scidb::Coordinate;
using scidb::Coordinates;

using scidb::ConstChunk;

using scidb::ConstArrayIterator;
using scidb::ConstItemIterator;
using scidb::ConstChunkIterator;

using scidb::Aggregate;
using scidb::AggregatePtr;
using scidb::AggregateLibrary;

using scidb::TypeLibrary;
using scidb::Type;
using scidb::Value;

using scidb::RLEPayload;

/**
 *  Pointer for an Array (very common in SciDb).
 */
typedef boost::shared_ptr<Array> ArrayPtr;

/**
 * A vector of Values
 */
typedef std::vector<Value> ValueVector;

/**
 * ConstItemIterator pointer
 */
typedef boost::shared_ptr<ConstItemIterator> ConstItemIteratorPtr;

#endif /* SEARCHLIGHT_SCIDB_INC_H_ */
