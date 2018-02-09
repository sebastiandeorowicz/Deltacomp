// ***********************************************************************************************
// The file is part of DeltaComp compressor
// Author  : Sebastian Deorowicz
// Version : 1.0
// Date    : 2018-01-10
// LIcense : GNU GPL3
// ***********************************************************************************************

#include "delta_comp.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <chrono>
#include <tuple>
#include <set>
#include <iterator>

using namespace std::chrono;

// ***********************************************************************************************
// DeltaComp constructor
CDeltaComp::CDeltaComp(int _max_delta_level, bool _verbose, int _fixed_delta_level, int _ppmd_order)
{
	max_delta_level = _max_delta_level;
	verbose = _verbose;

	fixed_delta_level = _fixed_delta_level;
	if (_ppmd_order >= 2)
	{
		ppmd_order = _ppmd_order;
		fixed_ppmd_order = true;
	}
	else
	{
		ppmd_order = PPMD_ORDER;
		fixed_ppmd_order = false;
	}

	vv_delta.resize(max_delta_level+1);

	ppmd_out = nullptr;

	lut_ilog2abs = new int[LUT_ILOG2ABS_SIZE];

	for (int i = 0; i < LUT_ILOG2ABS_SIZE; ++i)
		lut_ilog2abs[i] = ilog2(i);
}

// ***********************************************************************************************
// DeltaComp destructor
CDeltaComp::~CDeltaComp()
{
	if (ppmd_out)
		delete[] ppmd_out;

	delete[] lut_ilog2abs;
}

// ***********************************************************************************************
// Store unsigned integer in vector 
void CDeltaComp::store_uint(vector<uint8_t> &v, uint64_t x, int no_bytes)
{
	for (int i = 0; i < no_bytes; ++i)
	{
		v.push_back(x & 0xff);
		x >>= 8;
	}
}

// ***********************************************************************************************
// Read unsigned integer from vector
void CDeltaComp::read_uint(vector<uint8_t> &v, uint64_t pos, uint64_t &x, int no_bytes)
{
	x = 0;

	for (int i = 0; i < no_bytes; ++i)
		x += ((uint64_t)v[pos + i]) << (8 * i);
}

// ***********************************************************************************************
// Compress series of integers 
bool CDeltaComp::Compress(vector<int64_t> &series, vector<uint8_t> &comp_series, int &_best_delta_level)
{
	// Calculation of necessary number of delta values
	if(verbose)
		cout << "Deltas calculation\n";

	high_resolution_clock::time_point t1, t2;
	
	calculate_deltas(series);

	// Look for the best delta
	if (verbose)
		cout << "Find best delta\n";
	find_best_delta_level();

	// PPMD compress the vector
	if (verbose)
		cout << "PPMD compress\n";
	ppmd_compress();

	// Store PPMD output in the output vector
	for (int i = 0; i < (int) ppmd_compressed_size; ++i)
		comp_series.push_back(ppmd_out[i]);

	// Store metadata - sizes, delta
	store_uint(comp_series, byte_code_size, sizeof(uint64_t));
	store_uint(comp_series, series.size() * sizeof(uint64_t), sizeof(uint64_t));
	store_uint(comp_series, best_delta_level, 1);

	_best_delta_level = best_delta_level;

	return true;
}

// ***********************************************************************************************
// Decompress vector
bool CDeltaComp::Decompress(vector<uint8_t> &comp_series, vector<int64_t> &series)
{
	comp_data = comp_series;

	orig_no_rows = (int) series.size();

	if (verbose)
		cout << "PPMD decompress\n";

	// Read metadata - sizes, delta
	read_uint(comp_series, comp_series.size() - 2 * sizeof(uint64_t) - 1, byte_code_size, sizeof(uint64_t));
	read_uint(comp_series, comp_series.size() - sizeof(uint64_t) - 1, ppmd_uncompressed_size, sizeof(uint64_t));
	
	uint64_t tmp;
	read_uint(comp_series, comp_series.size() - 1, tmp, 1);
	best_delta_level = (int) tmp;

	// PPMD decompress of the compressed vector
	ppmd_compressed_size = comp_series.size() - 2 * sizeof(uint64_t);
	ppmd_decompress();

	// Perform undelta
	undelta();

	// Store the uncompressed data in the output vector
	series.clear();
	series.reserve(d_deltas.size());
	for (auto x : d_deltas)
		series.push_back(x);

	return true;
}

// ***********************************************************************************************
// Return best delta level
int CDeltaComp::GetDeltaLevel()
{
	return best_delta_level;
}

// ***********************************************************************************************
// Calculate integer log_2
int CDeltaComp::ilog2(int64_t x)
{
	int r;

	if (x == 0)
		return 0;
	else if (x == 1)
		return 1;
	else if (x <= 3)
		return 2;

	x >>= 2;

	for (r = 2; x; ++r)
		x >>= 1;

	return r;
}

// ***********************************************************************************************
// Calculate integer log_2 of absolute value
int CDeltaComp::ilog2abs(int64_t x)
{
	if (x < 0)
		x = -x;

	if (x < LUT_ILOG2ABS_SIZE)
		return lut_ilog2abs[x];

	int r = LOG_LUT_ILOG2ABS_SIZE;

	x >>= LOG_LUT_ILOG2ABS_SIZE;

	for (; x; ++r)
		x >>= 1;

	return r;
}

// ***********************************************************************************************
// Calculate integer absolute value
int64_t CDeltaComp::iabs(int64_t x)
{
	if (x < 0)
		return -x;
	else
		return x;
}

// ***********************************************************************************************
// Calculate deltas for series of integers
void CDeltaComp::calculate_deltas(vector<int64_t> &series)
{
	auto series_len = series.size();

	vv_delta[0].assign(series.begin(), series.end());
	v_hist_log2.resize(max_delta_level + 1);
	v_hist_log2[0].resize(65, 0);

	if (verbose)
		cout << 0 << " ";

	// Calculate entropy of the original data
	for (int j = 0; j < (int) series_len; ++j)
		++v_hist_log2[0][ilog2abs(vv_delta[0][j])];

	// Calculate deltas
	for (int i = 1; i <= max_delta_level; ++i)
	{
		if (verbose)
			cout << i << " ";
		v_hist_log2[i].resize(65, 0);

		auto &prev_vv_delta = vv_delta[i-1];
		auto &cur_vv_delta = vv_delta[i];
		auto &cur_v_hist_log2 = v_hist_log2[i];

		cur_vv_delta.resize(series_len);
		cur_vv_delta[0] = vv_delta[i - 1].front();
		++cur_v_hist_log2[ilog2abs(cur_vv_delta.front())];

		for (int j = 1; j < (int) series_len; ++j)
		{
			int64_t d = prev_vv_delta[j] - prev_vv_delta[j - 1];
			cur_vv_delta[j] = d;
			++cur_v_hist_log2[ilog2abs(d)];
		}
	}

	if (verbose)
		cout << endl;
}

// ***********************************************************************************************
// Find the best level of delta
void CDeltaComp::find_best_delta_level() 
{
	vector<double> v_entropy(max_delta_level + 1, 0.0);
	double sum = (double) vv_delta[0].size();

	for (int i = 0; i <= max_delta_level; ++i)
	{
		double entr = 0.0;

		for (int j = 0; j < (int) v_hist_log2[i].size(); ++j)
		{
			double cnt = (double) v_hist_log2[i][j];
			if (cnt == 0.0)
				continue;

			entr += -cnt * log2(cnt / sum) + cnt * j;
		}

		v_entropy[i] = entr;
		if (verbose)
			cout << i << " : " << (int64_t) entr << endl;
	}

	// Find best level
	auto p = min_element(v_entropy.begin(), v_entropy.end());
	best_delta_level = (int) distance(v_entropy.begin(), p);

	// Just for test of influence of the delta level
	if (fixed_delta_level >= 0)
		best_delta_level = fixed_delta_level;

	if (verbose)
		cout << best_delta_level << endl;
}

// ***********************************************************************************************
// Perform undelta
void CDeltaComp::undelta()
{
	for (int i = best_delta_level - 1; i >= 0; --i)
		for (int j = 1; j < (int) d_deltas.size(); ++j)
			d_deltas[j] = d_deltas[j - 1] + d_deltas[j];
}

// ***********************************************************************************************
// PPMD compress of the buffer for single level of delta and single PPMD order
pair<unsigned char*, uint64_t> CDeltaComp::ppmd_compress_single_level_single_order(int delta_level, int order)
{
	PpmdEncoder ppmd_e;
	uint64_t max_in_size = vv_delta[delta_level].size() * sizeof(uint64_t);
	byte_code_size = 0;

	unsigned char* ppmd_in = new unsigned char[max_in_size];

	// ByteCode the best deltas
	for (auto x : vv_delta[delta_level])
	{
		auto x_abs = iabs(x);
		int is_neg = x < 0;

		if (x_abs < 120)
			ppmd_in[byte_code_size++] = (unsigned char)(x_abs * 2 - is_neg);		// 0 -> 0, -1 -> 1, +1 -> 2, -2 -> 3, +2 -> 4, ..., +119 -> 238
		else
		{
			int x_log2 = ilog2(x_abs);
			int no_bytes = x_log2 / 8 + 1;
			int seg = 2 * no_bytes - is_neg + 238;
			ppmd_in[byte_code_size++] = seg;

			for (int i = 0; i < no_bytes; ++i)
			{
				ppmd_in[byte_code_size++] = x_abs & 0xff;
				x_abs >>= 8;
			}
		}
	}

	// PPMD compress of the byte-coded result
	uint64_t out_buf_size = byte_code_size + 1024;
	unsigned char* ppmd_out = new unsigned char[out_buf_size];
	ppmd_compressed_size = out_buf_size;

	ppmd_e.Encode(ppmd_in, byte_code_size, ppmd_out, ppmd_compressed_size, order, PPMD_MEMORY);

	if (verbose)
		cout << "PPMD compress: DL: " << delta_level << "    order: " << order << "    PPMD: " << ppmd_compressed_size << endl;

	delete[] ppmd_in;

	return make_pair(ppmd_out, ppmd_compressed_size);
}

// ***********************************************************************************************
// PPMD compress of the buffer
void CDeltaComp::ppmd_compress()
{
	if (fixed_delta_level >= 0 && fixed_ppmd_order)
	{
		tie(ppmd_out, ppmd_compressed_size) = ppmd_compress_single_level_single_order(fixed_delta_level, ppmd_order);
		if (verbose)
			cout << "Delta level: " << best_delta_level << "   PPMD order: " << ppmd_order << "      PPMD size: " << ppmd_compressed_size << endl;
		return;
	}

	unsigned char *best_out = nullptr;
	unsigned char *cur_out = nullptr;
	uint64_t best_size = 1ull << 61;
	uint64_t cur_size = 1ull << 61;

	// Adjust delta level
	if (fixed_delta_level < 0)
	{
		int order = fixed_ppmd_order ? ppmd_order : PPMD_ORDER;
		set<int> tested_delta_levels;

		tie(best_out, best_size) = ppmd_compress_single_level_single_order(best_delta_level, order);
		tested_delta_levels.insert(best_delta_level);

		if (verbose)
			cout << "Delta level: " << best_delta_level << "   PPMD order: " << order << "      PPMD size: " << best_size << endl;

		while (true)
		{
			if (best_delta_level > 1 && !tested_delta_levels.count(best_delta_level - 1))
			{
				int cur_delta_level = best_delta_level - 1;
				tie(cur_out, cur_size) = ppmd_compress_single_level_single_order(cur_delta_level, order);
				tested_delta_levels.insert(cur_delta_level);

				if (verbose)
					cout << "Delta level: " << cur_delta_level << "   PPMD order: " << order << "      PPMD size: " << cur_size << endl;

				if (cur_size < best_size)
				{
					delete[] best_out;
					best_out = cur_out;
					best_size = cur_size;
					best_delta_level = cur_delta_level;
				}
				else
					delete[] cur_out;
			}
			else if (best_delta_level < max_delta_level && !tested_delta_levels.count(best_delta_level + 1))
			{
				int cur_delta_level = best_delta_level + 1;
				tie(cur_out, cur_size) = ppmd_compress_single_level_single_order(cur_delta_level, order);
				tested_delta_levels.insert(cur_delta_level);

				if (verbose)
					cout << "Delta level: " << cur_delta_level << "   PPMD order: " << order << "      PPMD size: " << cur_size << endl;

				if (cur_size < best_size)
				{
					delete[] best_out;
					best_out = cur_out;
					best_size = cur_size;
					best_delta_level = cur_delta_level;
				}
				else
					delete[] cur_out;
			}
			else
				break;
		}
	}
	else
	{
		int order = fixed_ppmd_order ? ppmd_order : PPMD_ORDER;

		tie(best_out, best_size) = ppmd_compress_single_level_single_order(fixed_delta_level, order);
	}

	// Adjust PPMD order
	if (!fixed_ppmd_order)
	{
		set<int> candidate_orders = { 4, 6, 8, 10, 12, 14, 16, 24, 32, 48, 64 };
		set<int> tested_orders;
		int best_order = ppmd_order;

		tested_orders.insert(best_order);

		while (true)
		{
			if (best_order != *candidate_orders.begin() && !tested_orders.count(*(--candidate_orders.find(best_order))))
			{
				int cur_order = *(--candidate_orders.find(best_order));
				tie(cur_out, cur_size) = ppmd_compress_single_level_single_order(best_delta_level, cur_order);
				tested_orders.insert(cur_order);

				if (verbose)
					cout << "Delta level: " << best_delta_level << "   PPMD order: " << cur_order << "      PPMD size: " << cur_size << endl;

				if (cur_size < best_size)
				{
					delete[] best_out;
					best_out = cur_out;
					best_size = cur_size;
					best_order = cur_order;
				}
				else
					delete[] cur_out;
			}
			else if (best_order != *candidate_orders.rbegin() && !tested_orders.count(*(++candidate_orders.find(best_order))))
			{
				int cur_order = *(++candidate_orders.find(best_order));
				tie(cur_out, cur_size) = ppmd_compress_single_level_single_order(best_delta_level, cur_order);
				tested_orders.insert(cur_order);

				if (verbose)
					cout << "Delta level: " << best_delta_level << "   PPMD order: " << cur_order << "      PPMD size: " << cur_size << endl;

				if (cur_size < best_size)
				{
					delete[] best_out;
					best_out = cur_out;
					best_size = cur_size;
					best_order = cur_order;
				}
				else
					delete[] cur_out;
			}
			else
				break;
		}
	}

	ppmd_out = best_out;
	ppmd_compressed_size = best_size;
}

// ***********************************************************************************************
// PPMD decompression
void CDeltaComp::ppmd_decompress()
{
	PpmdDecoder ppmd_d;

	ppmd_out = new unsigned char[byte_code_size];
	unsigned char* ppmd_in = new unsigned char[ppmd_compressed_size];

	for (int i = 0; i < (int) ppmd_compressed_size; ++i)
		ppmd_in[i] = comp_data[i];

	ppmd_uncompressed_size = byte_code_size;
	ppmd_d.Decode(ppmd_in, ppmd_compressed_size, ppmd_out, ppmd_uncompressed_size, PPMD_MEMORY);

	if (verbose)
		cout << "Byte code size: " << ppmd_uncompressed_size << endl;

	// Decode byte code
	d_deltas.clear();

	for (int i = 0; i < (int) byte_code_size;)
	{
		int b = ppmd_out[i++];

		if (b <= 238)
		{
			if (b & 1)
				d_deltas.push_back(-(b + 1) / 2);
			else
				d_deltas.push_back(b / 2);			
		}
		else
		{
			int is_neg = b & 1;
			int no_bytes = (b + 1) / 2 - 119;
			uint64_t x = 0;

			for (int j = 0; j < no_bytes; ++j)
				x += ((uint64_t)ppmd_out[i++]) << (8 * j);

			if (is_neg)
				d_deltas.push_back(-(int64_t) x);
			else
				d_deltas.push_back(x);
		}
	}

	delete[] ppmd_in;
}

// EOF

