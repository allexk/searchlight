/* Copyright 2015, Brown University, Providence, RI.
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
 * @file dft_sampler.cpp
 *
 * DFT sampler for one-dimensional array.
 *
 * The input is assumed to be in CSV format with two columns: index and
 * value. Indices start from 0. If the value for an index is missing, it
 * is assumed to be empty (in SciDB sense). The output is in CSV format:
 * index, coordinate, low, high; where:
 * 	-- index - MBR index
 * 	-- coordinate -- DFT component number (0 through param)
 * 	-- low - low MBR coordinate
 * 	-- high - high MBR coordinate
 *
 * Params:
 * 	-- omega (-w 100) - subsequence size
 * 	-- MBR (-m 100) - MBR size (number of subsequences)
 * 	-- DFT (-d 3) - number of DFT components (complex numbers)
 * 	-- filename - file to process
 *
 * The output is given to the standard output.
 *
 * @author Alexander Kalinin
 */

#include <assert.h>
#include <unistd.h>
#include <fftw3.h>

#include <vector>
#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

/**
 * Coordinates type.
 */
typedef std::vector<double> Coordinates;

/**
 * MBR class for a bunch of points.
 */
class MBR {
public:
	/**
	 * Construct a new MBR.
	 *
	 * @param dims dimensionality
	 * @param mbr_size MBR size (in the number of sequences)
	 */
	MBR(size_t dims, size_t mbr_size) :
			low_(dims),
			high_(dims),
			mbr_size_(mbr_size) {
		Reset(0);
	}

	/**
	 * Resets MBR to cover the specified coordinate.
	 *
	 * @param coord coordinate to reset to
	 */
	void Reset(size_t coord) {
		start_pos_ = coord - coord % mbr_size_;
		end_pos_ = start_pos_ + mbr_size_ - 1;
		id_ = coord / mbr_size_;
		valid_ = false;
	}

	/**
	 * Check if the MBR is valid.
	 *
	 * @return true, if the MBR is valid
	 */
	bool Valid() const {
		return valid_;
	}

	/**
	 * Check if coordinate belongs to an MBR sequence.
	 *
	 * @param pos coordinate position
	 * @return true, if the coordinate belongs to the MBR
	 */
	bool CheckCoordinate(size_t pos) const {
		return pos >= start_pos_ && pos <= end_pos_;
	}

	/**
	 * Add a new point to the MBR.
	 *
	 * @param point point to add
	 */
	void AddPoint(const Coordinates &point) {
		assert(point.size() == low_.size());
		if (!valid_) {
			low_ = point;
			high_ = point;
			valid_ = true;
		} else {
			for (size_t i = 0; i < point.size(); ++i) {
				if (low_[i] > point[i]) {
					low_[i] = point[i];
				} else if (high_[i] < point[i]) {
					high_[i] = point[i];
				}
			}
		}
	}

	/**
	 * Output this MBR to the stream in CSV format.
	 *
	 * The format is:
	 * 	id,coord,low,high -- for every coordinate
	 *
	 * If the MBR is not valid, nothing is output.
	 *
	 * @param s stream to output to
	 */
	void OutputCSV(std::ostream &s) const {
		if (!valid_) return;

		for (size_t i = 0; i < low_.size(); ++i) {
			s << id_ << ',' << i << ',' << low_[i] << ',' << high_[i] << '\n';
		}
	}

private:
	int start_pos_ = 0; // Start sequence position
	int end_pos_ = 0; // End sequence position (inclusive)
	int id_ = 0; // ID of the MBR
	Coordinates low_, high_; // Low/high coordinates
	const size_t mbr_size_; // MBR size
	bool valid_ = false; // true if it's valid
};

/**
 * A very simple CSV reader.
 */
class CSVReader {
public:
	/**
	 * Create a new CSV reader from file.
	 *
	 * @param filename file name
	 */
	CSVReader(const std::string &filename) :
		inf_(filename) {}

	/**
	 * Parse next line.
	 *
	 * The assumed format is:
	 * 	coordinate,value
	 *
	 * @return true, if the value was read; false, otherwise
	 */
	bool Next() {
		std::string line;
		if (std::getline(inf_, line)) {
		    typedef boost::tokenizer<boost::char_separator<char>> tokenizer_t;
		    boost::char_separator<char> sep{","};
		    tokenizer_t tokenizer{line, sep};

		    tokenizer_t::const_iterator cit = tokenizer.begin();
		    if (cit == tokenizer.end()) {
		    	std::ostringstream err;
		    	err << "Cannot parse coordinate: '" << line << "'\n";
		    	throw std::runtime_error(err.str());
		    }
		    last_pos_ = boost::lexical_cast<size_t>(cit->c_str());
		    ++cit;
		    if (cit == tokenizer.end()) {
		    	std::ostringstream err;
		    	err << "Cannot parse value: '" << line << "'\n";
		    	throw std::runtime_error(err.str());
		    }
		    last_val_ = boost::lexical_cast<double>(cit->c_str());
		    ++cit;
		    if (cit != tokenizer.end()) {
		    	std::ostringstream err;
		    	err << "Unknown component in CSV: '" << line << "'\n";
		    	throw std::runtime_error(err.str());
		    }
		    return true;
		}
		return false;
	}

	/**
	 * Check if the reader is in EOF state.
	 *
	 * @return true, if the reader reached EOF.
	 */
	bool Eof() const {
		return inf_.eof();
	}

	/**
	 * Return last parsed coordinate.
	 *
	 * @return last parsed coordinate
	 */
	size_t Coordinate() const {
		return last_pos_;
	}

	/**
	 * Return last parsed value.
	 *
	 * @return last parsed value
	 */
	double Value() const {
		return last_val_;
	}

private:
	std::ifstream inf_;
	size_t last_pos_ = 0;
	double last_val_ = 0;
};

void SynchronizeMBR(MBR &mbr, size_t new_coordinate, std::ostream &s) {
	if (mbr.CheckCoordinate(new_coordinate)) {
		return;
	}
	// Either not valid or finished the current one
	if (mbr.Valid()) {
		// Finished a valid MBR
		mbr.OutputCSV(s);
	}
	mbr.Reset(new_coordinate);
}

void CheckAndHandleEOF(const CSVReader &reader, const MBR &mbr,
		std::ostream &s) {
	if (reader.Eof()) {
		mbr.OutputCSV(s);
		exit(0);
	}
}

int main(int argc, char **argv) {
	size_t w = 100;
	size_t m = 100;
	size_t d = 3;

	// Parse parameters
	int opt;
	while ((opt = getopt(argc, argv, "w:m:d:")) != -1) {
		int param;
		switch (opt) {
		case 'w':
			param = atoi(optarg);
			if (param <= 0) {
				std::cerr << "Omega must be positive!\n";
				return -1;
			}
			w = param;
			break;
		case 'm':
			param = atoi(optarg);
			if (param <= 0) {
				std::cerr << "MBR size must be positive!\n";
				return -1;
			}
			m = param;
			break;
		case 'd':
			param = atoi(optarg);
			if (param <= 0) {
				std::cerr << "DFT number size must be positive!\n";
				return -1;
			}
			d = param;
			break;
		default:
			std::cerr << "Unknown option: -" << opt;
			return -1;
		}
	}
	if (optind >= argc) {
		std::cerr << "Input filename is not specified\n";
		return -1;
	}
	const std::string in_file{argv[optind++]};

	// Open file and read it line by line
	CSVReader reader{in_file};
	std::ostream &outf = std::cout;
	MBR mbr(d * 2, m);
	// FFTW stuff
	double *in_fftw = (double *)fftw_malloc(sizeof(double) * w);
	fftw_complex *out_fftw = (fftw_complex *)fftw_malloc(
			sizeof(fftw_complex) * (w / 2 + 1));
	fftw_plan p = fftw_plan_dft_r2c_1d(w, in_fftw, out_fftw, FFTW_PATIENT);
	const double sqrt_N = std::sqrt(w);
	// Start parsing
	bool new_seq = true;
	size_t expected_coord = 0;
	size_t first_seq_coord = 0;
	size_t i = 0;
	while (true) {
		// Now we need to fill in the FFT buffer
		if (!new_seq) {
			first_seq_coord++;
			memmove(in_fftw, in_fftw + 1, sizeof(double) * (w - 1));
			i = w - 1;
		}
		bool empty_found = false;
		while (i < w) {
			if (!reader.Next()) {
				CheckAndHandleEOF(reader, mbr, outf);
			}
			size_t coord = reader.Coordinate();
			if (coord != expected_coord) {
				// Empty element
				expected_coord = coord + 1;
				in_fftw[0] = reader.Value();
				i = 1;
				first_seq_coord = coord;
				empty_found = true;
				break;
			} else {
				if (new_seq && i == 0) {
					first_seq_coord = coord;
				}
				in_fftw[i] = reader.Value();
				expected_coord++;
			}
			++i;
		}
		new_seq = empty_found;
		if (empty_found) continue;
		i = 0;
		// Can compute FFTW
		fftw_execute(p);
		Coordinates point;
		point.reserve(2 * d);
		for (size_t j = 0; j < d; ++j) {
			point.push_back(out_fftw[j][0] / sqrt_N);
			point.push_back(out_fftw[j][1] / sqrt_N);
		}
		// Add to the MBR
		SynchronizeMBR(mbr, first_seq_coord, outf);
		mbr.AddPoint(point);
	}

	// Cleanup
	fftw_destroy_plan(p);
	fftw_free(in_fftw);
	fftw_free(out_fftw);
	return 0;
}
