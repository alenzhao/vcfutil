#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include "algo.hpp"

std::string progname = "deep";


namespace vcf {

	using namespace std;
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	namespace ph = boost::phoenix;

	//FIXME Reader assumes all data lies on one chromosome

	class Reader {

	public:

		Reader(int window=0, bool fixed_win=false, int capacity=0)
			: fixed_window(fixed_win)
		{
			if (capacity > 0) {
				chromosomes.reserve(capacity);
				positions.reserve(capacity);
				values.reserve(capacity);
			}
			if (window < 0) {
				throw invalid_argument("Invalid window size");
			}
			if (window > 0) {
				win_chromosomes.reserve(window);
				win_positions.reserve(window);
				win_values.reserve(window);
			}
			this->window = window;
		}

		~Reader() {}

		void read(const string& filename) {
			ifstream file;
			file.open(filename.c_str(), ios::in);
			if (!file.is_open()) throw runtime_error("Failed to open input file.");
			read(file);
			file.close();
		}

		void read(istream& file) {
			string line;

			const char delim = '\t';
		
			// read in data
			size_t w = 0;
			// N.B. w has different meanings, contingent on fixed_window
			while (true) {
				getline(file, line);
				if (file.eof()) break;

				// ignore comment line
				if (line[0] == '#') continue;

				size_t parse_start = 0, parse_end;

				// first two elements are coordinates
				parse_end = find_nth_of(line, delim, 2);
				parse_coordinates(line.begin(), line.begin() + parse_end);

				// 8th element contains the information
				// start the end of the 7th element
				parse_start = find_nth_of(line, delim, 5, parse_end+1);
				parse_end = find_nth_of(line, delim, 1, parse_start+1);
				parse_depth(line.begin() + parse_start, line.begin() + parse_end);

				if (fixed_window) {
					// window size is fixed based on genomic coordinate: increment by the number of bp covered
					// N.B. one extra value from the next window will be summarized in the current window
					// w is index of the current window
					// check if the current position is greater than the limit of the current window
					if (win_positions.back() >= (w+1) * window) {
						// newly added data spills into the next window: summarize all elements in previous window

						// store the new element
						unsigned long position = win_positions.back();
						int chromosome = win_chromosomes.back();
						int value = win_values.back();

						// pop new element off the vector, temporarily
						win_positions.pop_back();
						win_chromosomes.pop_back();
						win_values.pop_back();

						summarize_data_fixed_window(chromosome, w*window + window/2);
						++w;

						// fast forward to the window containing the new element, summarizing data as necessary
						// each window spans [start_w, end_w)
						while( position >= (w+1) * window ) {
							summarize_data_fixed_window(chromosome, w*window + window/2);
							++w;
						}

						// push new element back onto the vector
						win_positions.push_back(position);
						win_chromosomes.push_back(chromosome);
						win_values.push_back(value);
					}
				} else {
					// w is a counter that represents number of elements accumulated in the current window
					// window size is based on cardinality of elements: increment counter by one element
					if (++w >= window) {
						// completed reading a window of data: summarize data
						summarize_data();
						w = 0;
					}
				}

			}

			// summarize the last partial window, if any
			if (fixed_window) {
				if (win_chromosomes.size() > 0) {
					summarize_data_fixed_window(win_chromosomes.front(), w*window + window/2);
				}
			} else {
				summarize_data();
			}

			// sanity check
			size = chromosomes.size();
			if (size != positions.size() || size != values.size()) {
				throw runtime_error("Input file is corrupt: array sizes are different");
			}

			// shrink as necessary
			chromosomes.resize(size);
			positions.resize(size);
			values.resize(size);
		}

		void print() {
			write(std::cout);
		}

		void write(const string& filename) {
			ofstream file;
			file.open(filename.c_str(), ios::out);
			if (!file.is_open()) throw runtime_error("Failed to open output file.");
			write(file);
			file.close();
		}

		void write(ostream& file) {
			const char delim = '\t';
			for (size_t i = 0; i < chromosomes.size(); ++i) {
				file << chromosomes[i] << delim << positions[i] << delim << values[i] << endl;
			}
		}

	private:

		template <typename iterator>
		bool parse_coordinates(iterator first, iterator last) {
			using qi::long_;
			using qi::_1;

			bool r = qi::phrase_parse(first, last,
				(
					long_[ph::push_back(ph::ref(win_chromosomes), _1)]
						>> long_[ph::push_back(ph::ref(win_positions), _1)]
				),
				ascii::space
			);

			if (first != last)  // partial match
				return false;
			return r;
		}

		template <typename iterator>
		bool parse_depth(iterator first, iterator last) {
			using qi::long_;
			using qi::_1;

			bool r = qi::phrase_parse(first, last,
				(
					"DP=" >> long_[ph::push_back(ph::ref(win_values), _1)]
				),
				ascii::space
			);
			if (first != last) // partial match
				return false;
			return r;
		}

		// find the nth instance of character c in string str
		size_t find_nth_of(const string& str, char c, size_t n = 1, size_t pos = 0) {
			size_t k = pos-1;
			for (size_t i = 0; i < n; ++i) {
				k = str.find_first_of(c, k+1);
			}
			return k;
		}

		// summarize data for current fixed window centered at given coordinates
		void summarize_data_fixed_window(int chromosome, unsigned long position) {
			if (win_values.size() == 0) {
				// window vectors is empty
				// add place-holder
				chromosomes.push_back( chromosome );
				positions.push_back( position );
				values.push_back( 0 );
			} else {
				// use the median to represent the data at the midpoint of the window	
				int median = algo::median( win_values.data(), win_values.size() );

				// push summarized data onto storing vectors
				chromosomes.push_back( chromosome );
				positions.push_back( position );
				values.push_back( median );

				// clear window vectors
				win_chromosomes.clear();
				win_positions.clear();
				win_values.clear();
			}
		}

		// summarize data and push onto storage vectors
		void summarize_data() {
			// window vector is empty: nothing to do
			if (win_values.size() == 0) return;

			// use the median to represent the data at the midpoint of the window	
			int median = algo::median( win_values.data(), win_values.size() );
			size_t middle = win_values.size() / 2;

			// push summarized data onto storing vectors
			chromosomes.push_back( win_chromosomes[middle] );
			positions.push_back( win_positions[middle] );
			values.push_back( median );

			// clear window vectors
			win_chromosomes.clear();
			win_positions.clear();
			win_values.clear();
		}


		// data within a window (temporary)
		vector<int> win_chromosomes;
		vector<unsigned long> win_positions;
		vector<int> win_values;

		// stored data (possibly downsampled)
		vector<int> chromosomes;
		vector<unsigned long> positions;
		vector<int> values;

		size_t size, window;

		bool fixed_window;

	};

}


int main(int argc, char *argv[]) {

	using namespace std;
	
	po::variables_map vm;
	po::options_description opts;
	po::positional_options_description popts;

	string input_fname, output_fname;
	int window;
	size_t length;
	bool fixed_window;

	// setup command line options

	opts.add_options()
		("help,h", "print help message")
		("input,i", po::value<string>(), "input file")
		("output,o", po::value<string>(&output_fname), "output file")
		("window,w", po::value<int>(&window)->default_value(0), "summarization window size")
		("fixed,f", po::value<bool>(&fixed_window)->default_value(true), "windows are fixed to genomic coordinate")
		("length,l", po::value<size_t>(&length)->default_value(1e6), "expected length of the values array")
	;
	popts.add("input", 1);

	po::store(
		po::command_line_parser(argc, argv)
			.options(opts)
			.positional(popts)
			.run(),
		vm
	);
	po::notify(vm);

	// process command line options

	if (vm.count("help")) {
		cout << "usage: " << progname << " [options] <input file> <output file>" << endl;
		cout << opts << endl;
		return 0;
	}

	if (vm.count("input")) {
		input_fname = vm["input"].as<string>();
	} else {
		throw invalid_argument("No input file specified.");
	}


	// read and write data

	vcf::Reader reader(window, fixed_window, length);
	reader.read(input_fname);

	if (output_fname.size() > 0) {
		reader.write(output_fname);
	} else {
		reader.print();
	}

	return 0;

}

