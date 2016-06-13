/* Copyright 2016, Brown University, Providence, RI.
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
 * @file paa_sampler.cpp
 *
 * PAA sampler for waveform data.
 *
 * @author Alexander Kalinin
 */


#include <gflags/gflags.h>
#include <boost/tokenizer.hpp>
#include "SciDBAPI.h"

#include <thread>
#include <exception>
#include <fstream>

// Command-line flags
DEFINE_int32(low_id, 0, "Low ID for waveform records");
DEFINE_int32(high_id, 99, "Low ID for waveform records");
DEFINE_int64(low_tick, 0, "Low tick for waveform records");
DEFINE_int64(high_tick, 50000000, "High tick for waveform records");
DEFINE_int32(seq_size, 128, "Sequence size for waveform subsequences");
DEFINE_int32(coeff, 8, "Number of PAA segments (coefficients)");
DEFINE_string(mbr, "1000,200,100", "MBR sizes for transformed PAAs");

using DataVector = std::vector<double>;

/**
 *  SciDB waveform reader.
 */
class SciDBReader {
public:
    /**
     * Construct a new reader.
     *
     * @param array array to read
     * @param attr attribute to read
     */
    SciDBReader(const std::string &array, const std::string &attr) :
                scidb_(scidb::getSciDB()),
                array_(array),
                attr_(attr) {
        std::cout << "Connecting to scidb...\n";
        conn_= scidb_.connect();
    }

    /**
     * Destructor.
     */
    ~SciDBReader() {
        std::cout << "Disconnecting from scidb...\n";
        scidb_.disconnect(conn_);
    }

    void ReadID(int id, DataVector &data) const {
        std::cout << "Reading waveform data for id=" << id << "...\n";
        // Length
        const size_t len = FLAGS_high_tick - FLAGS_low_tick + 1;
        std::cout << "Resizing data for " << len << " elements\n";
        data.resize(len);
        // Query to run (AFL)
        std::ostringstream query;
        query << "between(" << array_ << ", " << id << ", " << FLAGS_low_tick <<
                ", " << id << ", " << FLAGS_high_tick << ")";
        // Query SciDB
        std::cout << "Querying SciDB with: '" << query.str() << "'...\n";
        scidb::QueryResult query_res;
        scidb_.executeQuery(query.str(), true, query_res, conn_);
        // Find the attribute
        const auto &res_desc = query_res.array->getArrayDesc();
        const auto &attrs = res_desc.getAttributes(true);
        scidb::AttributeID attr_id = scidb::INVALID_ATTRIBUTE_ID;
        scidb::TypeId attr_type;
        for (const auto &desc: attrs) {
            if (desc.getName() == attr_) {
                attr_id = desc.getId();
                attr_type = desc.getType();
                break;
            }
        }
        if (attr_id == scidb::INVALID_ATTRIBUTE_ID) {
            throw std::runtime_error("Cannot find attribute " + attr_);
        }
        // Retrieve the data
        auto iter = query_res.array->getItemIterator(attr_id);
        while (!iter->end()) {
            const auto &pos = iter->getPosition();
            // Assume (id, tick)
            assert(pos[0] == id && pos[1] >= 0 &&
                   pos[1] - FLAGS_low_tick < data.size());
            data[pos[1]] = scidb::ValueToDouble(attr_type, iter->getItem());
            ++(*iter);
        }
        std::cout << "Finished reading waveform\n";
        scidb_.completeQuery(query_res.queryID, conn_);
    }

private:
    // SciDB client instance
    const scidb::SciDB &scidb_;
    // SciDB connection
    void *conn_;
    // Array to read
    const std::string array_;
    // Attribute to read
    const std::string attr_;
};

/**
 * Samples waveforms as PAA and outputs MBRs.
 */
class Sampler {
public:
    /**
     * Create a new sampler.
     *
     * @param array array name
     * @param attr attribute name
     */
    Sampler(const std::string &array, const std::string &attr) :
            reader_(array, attr) {
        // Parse MBR params
        using TokenSeparator = boost::char_separator<char>;
        using Tokenizer = boost::tokenizer<TokenSeparator>;
        TokenSeparator sep(",|"); // size_1xsize_2x...xsize_n
        Tokenizer tokenizer(FLAGS_mbr, sep);
        for (auto cit = tokenizer.begin(); cit != tokenizer.end(); ++cit) {
            const int mbr = std::stoi(*cit);
            std::ostringstream mbr_file;
            mbr_file << array << '_' << attr << '_' << FLAGS_seq_size <<
                    '_' << FLAGS_coeff << '_' << mbr << ".csv";
            mbr_files_.emplace_back(mbr,
                    File(new std::ofstream(mbr_file.str())));
            std::cout << "Created file " << mbr_file.str() <<
                    " for MBR size " << mbr << '\n';
        }
        if (mbr_files_.empty()) {
            throw std::invalid_argument("MBR parameter is wrong. "
                    "Cannot parse any sizes: " + FLAGS_mbr);
        }
    }

    /**
     * Process all ids.
     */
    void Process() {
        for (int32_t id = FLAGS_low_id; id <= FLAGS_high_id; ++id) {
            std::cout << "Processing id=" << id << "...\n";
            reader_.ReadID(id, data_);
            ProcessData(id);
        }
    }

private:
    // Types
    using File = std::shared_ptr<std::ostream>;
    using OutPair = std::pair<size_t, File>;

    // MBR processor
    class MBRProcessor {
    public:
        // Construct a new processor.
        MBRProcessor(size_t mbr_size, int id, DataVector &data,
                     const File &file) :
                mbr_size_(mbr_size),
                id_(id),
                data_(data),
                file_(file) {}

        // Execute the processor
        void Execute() {
            thr_ = std::thread(std::ref(*this));
        }

        // Call operator.
        void operator()() const {
            std::cout << "Processing data for MBR " << mbr_size_ << '\n';
            Process();
        }

        // Wait for processor
        void Wait() {
            if (thr_.joinable()) {
                // Blocking call
                thr_.join();
            }
        }

    private:
        // Process the data
        void Process() const {
            const size_t paa_size = FLAGS_seq_size / FLAGS_coeff;
            size_t mbr_num = 0;
            size_t mbr_count = 0;
            DataVector low(FLAGS_coeff), high(FLAGS_coeff), curr(FLAGS_coeff);
            for (size_t i = 0; i + FLAGS_seq_size - 1 < data_.size(); ++i) {
                for (size_t j = 0; j < FLAGS_coeff; ++j) {
                    curr[j] = data_[i + j * paa_size];
                }
                if (mbr_count == 0) {
                    low = high = curr;
                } else {
                    for (size_t j = 0; j < FLAGS_coeff; ++j) {
                        low[j] = std::min(low[j], curr[j]);
                        high[j] = std::max(high[j], curr[j]);
                    }
                }
                ++mbr_count;
                if (mbr_count == mbr_size_) {
                    OutputMBR(low, high, mbr_num);
                    ++mbr_num;
                    mbr_count = 0;
                }
            }
            // Last incomplete MBR (if started)
            if (mbr_count > 0) {
                OutputMBR(low, high, mbr_num);
            }
        }

        // Output MBR to file
        void OutputMBR(const DataVector &low, const DataVector &high,
                       size_t mbr_num) const {
            // Format: id, mbr_num, coeff, low_val, high_val
            std::ostream &f = *file_.get();
            for (size_t i = 0; i < low.size(); ++i) {
                f << id_ << ',' << mbr_num << ',' << i << ',' << low[i] <<
                        ',' << high[i] << '\n';
            }
        }

        const size_t mbr_size_;
        const int id_;
        DataVector &data_;
        File file_;
        std::thread thr_;
    };

    // Processes the current data array
    void ProcessData(int id) {
        assert(FLAGS_seq_size % FLAGS_coeff == 0);
        const size_t window_size = FLAGS_seq_size / FLAGS_coeff;
        assert(data_.size() >= window_size);
        // First we need to go with a sliding window
        double win_sum = 0.0;
        for (size_t i = 0; i < window_size - 1; ++i) {
            win_sum += data_[i];
        }
        for (size_t first = 0, last = window_size - 1; last < data_.size();
                ++first, ++last) {
            win_sum += data_[last];
            const double curr_avg = win_sum / window_size;
            win_sum -= data_[first];
            data_[first] = curr_avg;
        }
        // Then, compute and output PAA coefficients (mutli-thread)
        std::vector<MBRProcessor> procs;
        for (const auto &op: mbr_files_) {
            procs.emplace_back(op.first, id, data_, op.second);
        }
        // Execute
        std::for_each(procs.begin(), procs.end(),
                      [](MBRProcessor &p){p.Execute();});
        // Wait
        std::for_each(procs.begin(), procs.end(),
                      [](MBRProcessor &p){p.Wait();});
    }

    // Data reader
    SciDBReader reader_;
    // Data buffer
    DataVector data_;
    // Output files
    std::vector<OutPair> mbr_files_;
};

int main(int argc, char **argv) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    // Check some flags
    if (FLAGS_low_tick < 0) {
        throw std::invalid_argument("Low tick cannot be less than 0!");
    }
    if (FLAGS_seq_size % FLAGS_coeff != 0) {
        throw std::invalid_argument("Sequence size must be divisible by the "
                "number of coefficients!");
    }
    if (argc < 3) {
        throw std::invalid_argument("Array or attribute are not specified!");
    }
    // Process
    std::cout << "Processing array " << argv[1] << " and attribute " <<
            argv[2] << '\n';
    Sampler sampler(argv[1], argv[2]);
    sampler.Process();
    return 0;
}
