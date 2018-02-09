#pragma once

// ***********************************************************************************************
// The file is part of DeltaComp compressor
// Author  : Sebastian Deorowicz
// Version : 1.0
// Date    : 2018-01-10
// LIcense : GNU GPL3
// ***********************************************************************************************

#include <vector>
#include <cstdint>
#include "../libs/ppmd/PPMd.h"

using namespace std;

// ***********************************************************************************************
// CeltaComp class
class CDeltaComp
{
	const int PPMD_ORDER = 12;
	const int PPMD_MEMORY = 256*4;

	const int LOG_LUT_ILOG2ABS_SIZE = 10;
	const int LUT_ILOG2ABS_SIZE = 1 << LOG_LUT_ILOG2ABS_SIZE;

	vector<vector<int64_t>> vv_delta;
	int max_delta_level;
	int best_delta_level;
	int fixed_delta_level;
	int ppmd_order;
	bool fixed_ppmd_order;
	int *lut_ilog2abs;

	vector<vector<int64_t>> v_hist_log2;

	vector<uint8_t> comp_data;

	unsigned char *ppmd_out;
	uint64_t byte_code_size;
	uint64_t ppmd_uncompressed_size;
	uint64_t ppmd_compressed_size;

	bool verbose;

	// For decompression
	int orig_no_rows;
	vector<int64_t> d_deltas;

	inline int ilog2(int64_t x);
	inline int64_t iabs(int64_t x);
	inline int ilog2abs(int64_t x);

	void store_uint(vector<uint8_t> &v, uint64_t x, int no_bytes);
	void read_uint(vector<uint8_t> &v, uint64_t pos, uint64_t &x, int no_bytes);
		
	void calculate_deltas(vector<int64_t> &series);
	void find_best_delta_level();
	void undelta();

	pair<unsigned char*, uint64_t> ppmd_compress_single_level_single_order(int delta_level, int order);

	void ppmd_compress();
	void ppmd_decompress();

public:
	CDeltaComp(int _max_delta_level = 5, bool _verbose = false, int _fixed_delta_level = -1, int _ppmd_order = -1);
	~CDeltaComp();

	bool Compress(vector<int64_t> &series, vector<uint8_t> &comp_series, int &_best_delta_level);
	bool Decompress(vector<uint8_t> &comp_series, vector<int64_t> &series);
	int GetDeltaLevel();
};

// EOF
