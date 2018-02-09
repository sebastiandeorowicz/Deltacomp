// ***********************************************************************************************
// The file is part of DeltaComp compressor
// Author  : Sebastian Deorowicz
// Version : 1.0
// Date    : 2018-01-10
// LIcense : GNU GPL3
// ***********************************************************************************************

#include <CCfits/CCfits>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <cstdint>
#include <random>
#include <chrono>

#include "delta_comp.h"

//#define LOG_VALUES

#define DEVELOPER_MODE

using namespace CCfits;
using namespace std;
using namespace std::chrono;

#define NORM(x, mn, mx)	((x) < (mn) ? (mn) : (x) > (mx) ? (mx) : (x))

enum class in_format_t { double_t, int_t, uint_t, byte_t };
enum class processing_t {FITS_compress, FITS_decompress, text_compress, text_decompress};

// Parameters
double max_error = 0;
int max_delta_level = 10;
int best_delta_level = -1;
int fixed_delta_level = -1;
int ppmd_order = -1;
const int delta_level_min = 1;
const int delta_level_max = 32;
string in_name;
string out_name;
string hdu_name;
string column_name;
in_format_t in_format;
in_format_t out_format;
processing_t processing;
string in_type;
bool verbose = false;

// Input data
vector<double> v_d_data;
vector<int> v_i_data;
vector<uint8_t> v_b_data;
int orig_no_rows;

// Internal data
vector<int64_t> cnv_data;

// Compressed data
vector<uint8_t> comp_data;

// Functions
void usage();
bool parse_params(int argc, char **argv);
bool FITS_read_raw_data();
bool FITS_read_compressed_data();
bool text_read_raw_data();
bool text_read_compressed_data();

bool convert_input_data();
bool convert_input_data_double();
bool convert_input_data_int();
bool convert_input_data_uint();
bool convert_output_data();
bool convert_output_data_double();
bool convert_output_data_int();
bool convert_output_data_uint();

bool compress();
bool decompress();
bool FITS_store_compressed_file();
bool FITS_store_uncompressed_file();
bool text_store_compressed_file();
bool text_store_uncompressed_file();

void FITS_compression();
void FITS_decompression();
void text_compression();
void text_decompression();

// ***********************************************************************************************
// Show usage
void usage()
{
	cerr << "Usage: deltacomp <mode> [Options] <in_name> <out_name> <hdu_name> <column_name> \n";
	cerr << "   mode        - one of: fc (FITS compress), fd (FITS decompress), tc (text compress), td (text decompress)\n";
	cerr << "   in_name     - input file name (in FITS format)\n";
	cerr << "   out_name    - output file name (in FITS format)\n";
	cerr << "   hdu_name    - hdu name in the input file (only for FITS files)\n";
	cerr << "   column_name - column to compress name (only for FITS files)\n";
	cerr << "Options:\n";
	cerr << "   -me|--max-error <value> - maximal value of the error (in lossy mode)\n";
	cerr << "   -mdl|--max-delta-level <value> - maximal value of delta level (default: " << max_delta_level << ")\n";
#ifdef DEVELOPER_MODE
	cerr << "   -fdl|--fixed-delta-level <value> - fixed value of delta level (default: " << fixed_delta_level << ")\n";
	cerr << "   -po|--ppmd-order <value> - PPMD order (default: " << ppmd_order << ")\n";
	cerr << "   -v - be verbose\n";
#endif
}

// ***********************************************************************************************
// Parse execution parameters
bool parse_params(int argc, char **argv)
{
	if (argc < 4)
	{
		cerr << "Too few parameters\n";
		usage();

		return false;
	}

	if (string(argv[1]) == "fc")
		processing = processing_t::FITS_compress;
	else if (string(argv[1]) == "fd")
		processing = processing_t::FITS_decompress;
	else if (string(argv[1]) == "tc")
		processing = processing_t::text_compress;
	else if (string(argv[1]) == "td")
		processing = processing_t::text_decompress;
	else
	{
		cerr << "Invalid mode: " << string(argv[1]) << endl;
		usage();
		return false;
	}

	int i = 2;

	while (i + 1 < argc)
	{
		if (string(argv[i]) == "-v")
		{
			verbose = true;
			i += 1;
		}
		if (string(argv[i]) == "-me" || string(argv[i]) == "--max-error")
		{
			max_error = atof(argv[i + 1]);
			i += 2;
		}
		else if (string(argv[i]) == "-mdl" || string(argv[i]) == "--max-delta-level")
		{
			max_delta_level = atoi(argv[i + 1]);
			max_delta_level = NORM(max_delta_level, delta_level_min, delta_level_max);
			i += 2;
		}
#ifdef DEVELOPER_MODE
		else if (string(argv[i]) == "-fdl" || string(argv[i]) == "--fixed-delta-level")
		{
			fixed_delta_level = atoi(argv[i + 1]);
			i += 2;
		}
		else if (string(argv[i]) == "-po" || string(argv[i]) == "--ppmd-order")
		{
			ppmd_order = atoi(argv[i + 1]);
			i += 2;
		}
#endif
		else
			break;
	}

	if (processing == processing_t::FITS_compress || processing == processing_t::FITS_decompress)
	{
		if (i + 4 != argc)
		{
			usage();
			return false;
		}

		in_name = string(argv[i++]);
		out_name = string(argv[i++]);
		hdu_name = string(argv[i++]);
		column_name = string(argv[i]);
	}
	else
	{
		if (i + 2 != argc)
		{
			usage();
			return false;
		}

		in_name = string(argv[i++]);
		out_name = string(argv[i++]);
	}

	return true;
}

// ***********************************************************************************************
// Read input FITS file
bool FITS_read_raw_data()
{
	vector<string> hdus = { hdu_name };

	FITS inFile(in_name, Read, hdus, false);

	Column &column = inFile.extension(hdus[0]).column(column_name);

	if(verbose)
		cout << "Rows: " << column.rows() << endl;

	string column_format = column.format();

	in_type = column_format;

	if (column_format == "D")
	{
		in_format = in_format_t::double_t;
		column.read(v_d_data, 1, column.rows());
	}
	else if (column_format == "I" || column_format == "K" || column_format == "J")
	{
		in_format = in_format_t::int_t;
		column.read(v_i_data, 1, column.rows());
	}
	else if (column_format == "U" || column_format == "V")
	{
		in_format = in_format_t::uint_t;
		column.read(v_i_data, 1, column.rows());
	}
	else if (column_format == "B")
	{
		in_format = in_format_t::byte_t;
		column.read(v_b_data, 1, column.rows());
	}
	else
	{
		cerr << "Incompatibile column format\n";
		return false;
	}

#ifdef LOG_VALUES
	FILE *log = fopen("in.log", "wb");
	for (auto x : v_d_data)
		fprintf(log, "%20.17f\n", x);
	fclose(log);
#endif

	return true;
}

// ***********************************************************************************************
// Read compressed FITS file
bool FITS_read_compressed_data()
{
	vector<string> hdus = { hdu_name };

	FITS inFile(in_name, Read, hdus, false);

	Column &column = inFile.extension(hdus[0]).column(column_name);

	if(verbose)
		cout << "Rows: " << column.rows() << endl;

	string column_format = column.format();

	if (column_format != "1B")
	{
		cerr << "Wrong input file\n";
		return false;
	}

	column.read(comp_data, 1, column.rows());

	inFile.extension(hdus[0]).readKey("DC_DL", max_delta_level);
	inFile.extension(hdus[0]).readKey("DC_O_NR", orig_no_rows);
	inFile.extension(hdus[0]).readKey("DC_O_T", in_type);
	inFile.extension(hdus[0]).readKey("DC_ME", max_error);

	return true;
}

// ***********************************************************************************************
// Read input text file (simple format: each line contains a single double)
bool text_read_raw_data()
{
	FILE *in = fopen(in_name.c_str(), "rb");
	char s[1024];

	if (!in)
		return false;

	setvbuf(in, nullptr, _IOFBF, 64 << 20);

	v_d_data.clear();

	while (!feof(in))
	{
		fgets(s, 1024, in);
		if (feof(in))
			break;

		double v = atof(s);

		v_d_data.push_back(v);
	}

	fclose(in);

	return true;
}

// ***********************************************************************************************
// Read compressed text file
bool text_read_compressed_data()
{
	FILE *in = fopen(in_name.c_str(), "rb");

	if (!in)
		return false;

	fseek(in, 0, SEEK_END);
	size_t size = ftell(in) - sizeof(double);
	fseek(in, 0, SEEK_SET);

	fread(&max_error, sizeof(double), 1, in);

	comp_data.resize(size);
	fread(comp_data.data(), 1, comp_data.size(), in);

	fclose(in);

	return true;
}

// ***********************************************************************************************
// Convert input data to integers - choose proper method depending on the input format
bool convert_input_data()
{
	if (in_format == in_format_t::double_t)
		return convert_input_data_double();
	else if (in_format == in_format_t::int_t)
		return convert_input_data_int();
	else if (in_format == in_format_t::uint_t)
		return convert_input_data_uint();
	else if (in_format == in_format_t::byte_t)
		return convert_input_data_int();
	
	return false;
}

// ***********************************************************************************************
// Convert data during decompression from integers to proper format
bool convert_output_data()
{
	if (in_type == "D")
		return convert_output_data_double();
	else if (in_type == "I" || in_type == "K" || in_type == "J")
		return convert_output_data_int();
	else if (in_type == "U" || in_type == "V")
		return convert_output_data_uint();
	else if (in_type == "B")
		return convert_output_data_int();

	return false;
}

// ***********************************************************************************************
// Convert imput doubles to integers (devide by max_error)
bool convert_input_data_double()
{
	if (max_error <= 0)
	{
		cerr << "Invalid max_error value (" << max_error << "). Changed to 1.0\n";
		max_error = 1;
	}

	cnv_data.resize(v_d_data.size());

	for (int i = 0; i < (int) v_d_data.size(); ++i)
		cnv_data[i] = (int64_t) round(v_d_data[i] / max_error);

	return true;
}

// ***********************************************************************************************
// Convert input intger data to integers (divide by max_error)
bool convert_input_data_int()
{
	if (max_error <= 0)
	{
		cerr << "Invalid max_error value (" << max_error << "). Changed to 1.0\n";
		max_error = 1;
	}

	cnv_data.resize(v_i_data.size());

	for (int i = 0; i < (int) v_i_data.size(); ++i)
		cnv_data[i] = (int64_t) (v_i_data[i] / max_error);

	return true;
}

// ***********************************************************************************************
// Convert input unsinged integer data to integers (divide by max_error)
bool convert_input_data_uint()
{
	if (max_error <= 0)
	{
		cerr << "Invalid max_error value (" << max_error << "). Changed to 1.0\n";
		max_error = 1;
	}

	cnv_data.resize(v_i_data.size());

	for (int i = 0; i < (int) v_i_data.size(); ++i)
		cnv_data[i] = (int64_t)(v_i_data[i] / max_error);

	return true;
}

// ***********************************************************************************************
// Convert output data (integers to doubles)
bool convert_output_data_double()
{
	v_d_data.resize(cnv_data.size());

	for (int i = 0; i < (int) v_d_data.size(); ++i)
		v_d_data[i] = (double) cnv_data[i] * max_error;

	return true;
}

// ***********************************************************************************************
// Convert output data (integers to integers)
bool convert_output_data_int()
{
	v_i_data.resize(cnv_data.size());

	for (int i = 0; i < (int) v_i_data.size(); ++i)
		v_i_data[i] = cnv_data[i] * max_error;

	return true;
}

// ***********************************************************************************************
// Convert output data (integers to unsigned integers)
bool convert_output_data_uint()
{
	v_i_data.resize(cnv_data.size());

	for (int i = 0; i < (int) v_i_data.size(); ++i)
		v_i_data[i] = cnv_data[i] * max_error;

	return true;
}

// ***********************************************************************************************
// Perform compression using DeltaComp algorithm
bool compress()
{
	CDeltaComp dc(max_delta_level, verbose, fixed_delta_level, ppmd_order);

	dc.Compress(cnv_data, comp_data, best_delta_level);

	return true;
}

// ***********************************************************************************************
// Perform decompression using DeltaComp algorithm
bool decompress()
{
	CDeltaComp dc(max_delta_level, verbose);

	cnv_data.resize(orig_no_rows);

	dc.Decompress(comp_data, cnv_data);

	return true;
}

// ***********************************************************************************************
// Store compressed data in FITS format
bool FITS_store_compressed_file()
{
	FITS outFile(out_name, Write);

	vector<string> colName(1);
	vector<string> colForm(1);
	vector<string> colUnit(1);
	
	colName[0] = column_name;
	colForm[0] = "1B";
	colUnit[0] = "";

	string hduName(hdu_name);

	Table *newTable = outFile.addTable(hduName, comp_data.size(), colName, colForm, colUnit);
	newTable->column(colName[0]).write(comp_data, 1);
	
	// Copy existing metadata
	vector<string> hdus = { hdu_name };
	FITS inFile(in_name, Read, hdus, false);
	ExtHDU& oldTable = inFile.extension(hdus[0]);
	newTable->copyAllKeys(&oldTable);

	// Add new metadata describing the new format
	newTable->addKey("AUTHOR", "Sebastian Deorowicz", "");
	newTable->addKey("COMPR", "DeltaComp", "Data stored in DeltaComp format");
	newTable->addKey("DC_DL", best_delta_level, "Delta level used for compression in DeltaComp");
	newTable->addKey("DC_O_NR", (int) cnv_data.size(), "Number of rows in uncompressed file");
	newTable->addKey("DC_O_T", in_type, "Format of uncompressed column");
	newTable->addKey("DC_ME", max_error, "Maximal error in lossy compression");

	return true;
}

// ***********************************************************************************************
// Store uncompressed FITS file
bool FITS_store_uncompressed_file()
{
	FITS outFile(out_name, Write);

	vector<string> colName(1);
	vector<string> colForm(1);
	vector<string> colUnit(1);

	string hduName(hdu_name);
	colName[0] = column_name;
	colForm[0] = in_type;
	colUnit[0] = "";

	Table *newTable = outFile.addTable(hduName, cnv_data.size(), colName, colForm, colUnit);
	Column &column = outFile.extension(hdu_name).column(column_name);

	// Copy existing metadata
	vector<string> hdus = { hdu_name };
	FITS inFile(in_name, Read, hdus, false);
	ExtHDU& oldTable = inFile.extension(hdus[0]);
	newTable->copyAllKeys(&oldTable);


#ifdef LOG_VALUES
	FILE *log = fopen("out.log", "wb");
	for (auto x : v_d_data)
		fprintf(log, "%20.17f\n", x);
	fclose(log);
#endif

	if (in_type == "D")
	{
		column.write(v_d_data, 1);
	}
	else if (in_type == "I" || in_type == "K" || in_type == "J")
	{
		column.write(v_i_data, 1);
	}
	else if (in_type == "U" || in_type == "V")
	{
		in_format = in_format_t::uint_t;
		column.write(v_i_data, 1);
	}
	else if (in_type == "B")
	{
		in_format = in_format_t::byte_t;
		column.write(v_b_data, 1);
	}
	else
	{
		cerr << "Incompatibile column format\n";
		return false;
	}

	return true;
}

// ***********************************************************************************************
// Store compressed text file
bool text_store_compressed_file()
{
	FILE *out = fopen(out_name.c_str(), "wb");
	
	if (!out)
	{
		cerr << "Cannot open " << out_name << " file\n";
		return false;
	}

	fwrite(&max_error, sizeof(double), 1, out);
	fwrite(comp_data.data(), 1, comp_data.size(), out);

	fclose(out);

	return true;
}

// ***********************************************************************************************
// Store uncompressed text file
bool text_store_uncompressed_file()
{
	FILE *out = fopen(out_name.c_str(), "wb");

	if (!out)
	{
		cerr << "Cannot open " << out_name << " file\n";
		return false;
	}

	setvbuf(out, nullptr, _IOFBF, 64 << 20);
	
	for (auto x : v_d_data)
		fprintf(out, "%.12f\n", x);

	fclose(out);

	return true;
}

// ***********************************************************************************************
// Perform FITS file compression
void FITS_compression()
{
	if (!FITS_read_raw_data())
	{
		cerr << "Data read failed\n";
		return;
	}

	convert_input_data();
	compress();
	try {
		FITS_store_compressed_file();
	}
	catch (...)
	{
		cerr << "Cannot store FITS file\n";
	}
}

// ***********************************************************************************************
// Perform FITS file decompression
void FITS_decompression()
{
	if (!FITS_read_compressed_data())
	{
		cerr << "Data read failed\n";
		return;
	}

	decompress();
	convert_output_data();

	try {
		FITS_store_uncompressed_file();
	}
	catch(...)
	{
		cerr << "Cannot store FITS file\n";
	}
}

// ***********************************************************************************************
// Perform text file compression
void text_compression()
{
	if (!text_read_raw_data())
	{
		cerr << "Data read failed\n";
		return;
	}

	in_format = in_format_t::double_t;
	if (max_error == 0.0)
		max_error = 1.0;

	convert_input_data();
	compress();

	text_store_compressed_file();
}

// ***********************************************************************************************
// Perform text file decompression
void text_decompression()
{
	if (!text_read_compressed_data())
	{
		cerr << "Data read failed\n";
		return;
	}

	decompress();
	in_type = "D";
	convert_output_data();

	text_store_uncompressed_file();
}

// ***********************************************************************************************
// Main function
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
		return 0;

	high_resolution_clock::time_point t1, t2;

	t1 = high_resolution_clock::now();

	if (processing == processing_t::FITS_compress)
		FITS_compression();
	else if (processing == processing_t::FITS_decompress)
		FITS_decompression();
	else if (processing == processing_t::text_compress)
		text_compression();
	else if (processing == processing_t::text_decompress)
		text_decompression();

	t2 = high_resolution_clock::now();
	if (verbose)
		cout << "Total  time: " << duration_cast<duration<double>>(t2 - t1).count() << endl;
	t1 = t2;
	
	return 0;
}

// EOF
