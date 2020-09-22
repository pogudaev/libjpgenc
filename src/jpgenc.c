/*
Copyright (c) 2020 Alexander Pogudaev
This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:
1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgement in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
    Alexander Pogudaev
    pogudaev@yandex.ru
*/

#include "jpgenc.h"
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <arpa/inet.h>
#include <stdbool.h>

#include "jpgenc_const.h"
#include "jpgenc_definition.h"

// ---------------- STATIC -----------------------
static uint8_t *jpgenc_huff_get_code_lengths(uint8_t huffsize[], uint8_t const *bits);
static uint16_t *jpgenc_huff_get_codes(uint16_t codes[], uint8_t *huffsize, int64_t count);
static void jpgenc_huff_get_extended(uint8_t *out_ehuffsize,
                                     uint16_t *out_ehuffcode,
                                     uint8_t const *huffval,
                                     uint8_t *huffsize,
                                     uint16_t *huffcode, int64_t count);
static void jpgenc_huff_expand(JPGEncObject *state);
static void jpgenc_write(JPGEncObject *state, const void *data, size_t num_bytes, size_t num_elements);
static void jpgenc_write_bits(JPGEncObject *state,
                              uint32_t *bitbuffer, uint32_t *location,
                              uint16_t num_bits, uint16_t bits);
static void jpgenc_write_DHT(JPGEncObject *state,
                             uint8_t const *matrix_len,
                             uint8_t const *matrix_val,
                             JPGEncHuffmanTableClass ht_class,
                             uint8_t id);
static void jpgenc_write_DQT(JPGEncObject *state, uint8_t id);
static void jpgenc_fdct (float *data);



// Returns:
//  out[1] : number of bits
//  out[0] : bits
static void jpgenc_calculate_variable_length_int(int value, uint16_t out[2])
{
	int abs_val = value;

	if ( value < 0 ) {
		abs_val = -abs_val;
		--value;
	}

	out[1] = 1;

	while ( abs_val >>= 1 ) {
		++out[1];
	}

	out[0] = (uint16_t)(value & ((1 << out[1]) - 1));
}


static void jpgenc_fdct_pass(float *data, int pass)
{
	float tmp10, tmp11, tmp12, tmp13;
	float tmp[8];

	float z1, z2, z3, z4, z5, z11, z13;
	float *dataptr;
	int ctr;

	/* Pass 1: process rows. */

	int k = (!pass) ? 1 : 8;

	dataptr = data;

	for ( ctr = 7; ctr >= 0; ctr-- ) {
		for (int i = 0; i < 4; i++) {
			tmp[i] = dataptr[k * i] + dataptr[k * (7 - i)];
		}

		for (int i = 4; i < 8; i++) {
			tmp[i] = dataptr[k * (7 - i)] - dataptr[k * i];
		}

		/* Even part */

		tmp10 = tmp[0] + tmp[3];    /* phase 2 */
		tmp11 = tmp[1] + tmp[2];
		tmp12 = tmp[1] - tmp[2];
		tmp13 = tmp[0] - tmp[3];

		dataptr[k * 0] = tmp10 + tmp11; /* phase 3 */
		dataptr[k * 4] = tmp10 - tmp11;

		z1 = (tmp12 + tmp13) * ((float) 0.707106781); /* c4 */
		dataptr[k * 2] = tmp13 + z1; /* phase 5 */
		dataptr[k * 6] = tmp13 - z1;

		/* Odd part */

		tmp10 = tmp[4] + tmp[5];    /* phase 2 */
		tmp11 = tmp[5] + tmp[6];
		tmp12 = tmp[6] + tmp[7];

		/* The rotator is modified from fig 4-8 to avoid extra negations. */
		z5 = (tmp10 - tmp12) * ((float) 0.382683433); /* c6 */
		z2 = ((float) 0.541196100) * tmp10 + z5; /* c2-c6 */
		z4 = ((float) 1.306562965) * tmp12 + z5; /* c2+c6 */
		z3 = tmp11 * ((float) 0.707106781); /* c4 */

		z11 = tmp[7] + z3;        /* phase 5 */
		z13 = tmp[7] - z3;

		dataptr[k * 5] = z13 + z2; /* phase 6 */
		dataptr[k * 3] = z13 - z2;
		dataptr[k * 1] = z11 + z4;
		dataptr[k * 7] = z11 - z4;

		dataptr += (pass) ? 1 : 8;          /* advance pointer to next column */
	}
}


static void jpgenc_fdct (float *data)
{
	jpgenc_fdct_pass(data, 0);
	jpgenc_fdct_pass(data, 1);
}

static void jpgenc_write_DQT(JPGEncObject *state, uint8_t id)
{
	uint16_t DQT = htons(DQT_MARKER);
	jpgenc_write(state, &DQT, sizeof(uint16_t), 1);
	uint16_t len = htons(2 /*len*/ + 1 /*id*/ + 64 /*matrix*/);
	jpgenc_write(state, &len, sizeof(uint16_t), 1);
	assert(id < 4);
	uint8_t precision_and_id = id;  // 0x0000 8 bits | 0x00id
	jpgenc_write(state, &precision_and_id, sizeof(uint8_t), 1);
	// Write matrix
	jpgenc_write(state, one_matrix, sizeof(one_matrix), 1);
}

// Returns all code sizes from the BITS specification (JPEG C.3)
static uint8_t *jpgenc_huff_get_code_lengths(uint8_t huffsize[/*256*/], uint8_t const *bits)
{
	int k = 0;

	for ( int i = 0; i < 16; /* nothing */) {
		int si = i + 1;

		for ( int j = 0; j < bits[i]; ++j ) {
			huffsize[k++] = (uint8_t)(si);
		}

		huffsize[k] = 0;
		i = si;
	}

	return huffsize;
}

// Fills out the prefixes for each code.
static uint16_t *jpgenc_huff_get_codes(uint16_t codes[], uint8_t *huffsize, int64_t count)
{
	uint16_t code = 0;
	int k = 0;
	uint8_t sz = huffsize[0];

	while (true) {
		do {
			assert(k < count);
			codes[k++] = code++;
		} while (huffsize[k] == sz);

		if (huffsize[k] == 0) {
			return codes;
		}

		do {
			code = (uint16_t)(code << 1);
			++sz;
		} while ( huffsize[k] != sz );
	}
}


static void jpgenc_huff_get_extended(uint8_t *out_ehuffsize,
                                     uint16_t *out_ehuffcode,
                                     uint8_t const *huffval,
                                     uint8_t *huffsize,
                                     uint16_t *huffcode, int64_t count)
{
	int k = 0;

	do {
		uint8_t val = huffval[k];
		out_ehuffcode[val] = huffcode[k];
		out_ehuffsize[val] = huffsize[k];
		k++;
	} while ( k < count );
}

static void jpgenc_huff_expand(JPGEncObject *state)
{
	assert(state);

	state->ht_bits[JPGENC_LUMA_DC]   = default_ht_luma_dc_len;
	state->ht_bits[JPGENC_LUMA_AC]   = default_ht_luma_ac_len;
	state->ht_bits[JPGENC_CHROMA_DC] = default_ht_chroma_dc_len;
	state->ht_bits[JPGENC_CHROMA_AC] = default_ht_chroma_ac_len;

	state->ht_vals[JPGENC_LUMA_DC]   = default_ht_luma_dc;
	state->ht_vals[JPGENC_LUMA_AC]   = default_ht_luma_ac;
	state->ht_vals[JPGENC_CHROMA_DC] = default_ht_chroma_dc;
	state->ht_vals[JPGENC_CHROMA_AC] = default_ht_chroma_ac;

	// How many codes in total for each of LUMA_(DC|AC) and CHROMA_(DC|AC)
	int32_t spec_tables_len[4] = { 0 };

	for ( int i = 0; i < 4; ++i ) {
		for ( int k = 0; k < 16; ++k ) {
			spec_tables_len[i] += state->ht_bits[i][k];
		}
	}

	// Fill out the extended tables..
	uint8_t huffsize[4][257];
	uint16_t huffcode[4][256];

	for ( int i = 0; i < 4; ++i ) {
		assert (256 >= spec_tables_len[i]);
		jpgenc_huff_get_code_lengths(huffsize[i], state->ht_bits[i]);
		jpgenc_huff_get_codes(huffcode[i], huffsize[i], spec_tables_len[i]);
	}

	for ( int i = 0; i < 4; ++i ) {
		int64_t count = spec_tables_len[i];
		jpgenc_huff_get_extended(state->ehuffsize[i],
		                         state->ehuffcode[i],
		                         state->ht_vals[i],
		                         &huffsize[i][0],
		                         &huffcode[i][0], count);
	}
}

static void jpgenc_write(JPGEncObject *state, const void *data, size_t num_bytes, size_t num_elements)
{
	size_t to_write = num_bytes * num_elements;

	// Cap to the buffer available size and copy memory.
	size_t capped_count = JPGENC_BUFFER_SIZE - 1 - state->output_buffer_count;

	if (capped_count > to_write) {
		capped_count = to_write;
	}

	memcpy(state->output_buffer + state->output_buffer_count, data, capped_count);
	state->output_buffer_count += capped_count;

	assert (state->output_buffer_count <= JPGENC_BUFFER_SIZE - 1);

	// Flush the buffer.
	if ( state->output_buffer_count == JPGENC_BUFFER_SIZE - 1 ) {
		state->write_context.func(state->write_context.context, state->output_buffer, state->output_buffer_count);
		state->output_buffer_count = 0;
	}

	// Recursively calling ourselves with the rest of the buffer.
	if (capped_count < to_write) {
		jpgenc_write(state, (const uint8_t *) data + capped_count, to_write - capped_count, 1);
	}
}


static void jpgenc_write_bits(JPGEncObject *state,
                              uint32_t *bitbuffer, uint32_t *location,
                              uint16_t num_bits, uint16_t bits)
{
	//   v-- location
	//  [                     ]   <-- bit buffer
	// 32                     0
	//
	// This call pushes to the bitbuffer and saves the location. Data is pushed
	// from most significant to less significant.
	// When we can write a full byte, we write a byte and shift.

	// Push the stack.
	uint32_t nloc = *location + num_bits;
	*bitbuffer |= (uint32_t)(bits << (32 - nloc));
	*location = nloc;

	while ( *location >= 8 ) {
		// Grab the most significant byte.
		uint8_t c = (uint8_t)((*bitbuffer) >> 24);
		// Write it to file.
		jpgenc_write(state, &c, 1, 1);

		if ( c == 0xff )  {
			// Special case: tell JPEG this is not a marker.
			char z = 0;
			jpgenc_write(state, &z, 1, 1);
		}

		// Pop the stack.
		*bitbuffer <<= 8;
		*location -= 8;
	}
}


static void jpgenc_write_DHT(JPGEncObject *state,
                             uint8_t const *matrix_len,
                             uint8_t const *matrix_val,
                             JPGEncHuffmanTableClass ht_class,
                             uint8_t id)
{
	int num_values = 0;

	for ( int i = 0; i < 16; ++i ) {
		num_values += matrix_len[i];
	}

	assert(num_values <= UINT16_MAX);

	uint16_t marker = htons(HAFF_TAB_MARKER);
	// 2(len) + 1(Tc|th) + 16 (num lengths) + ?? (num values)
	uint16_t len = htons(2 + 1 + 16 + (uint16_t)num_values);
	assert(id < 4);
	uint8_t tc_th = (uint8_t)((((uint8_t)ht_class) << 4) | id);

	jpgenc_write(state, &marker, sizeof(uint16_t), 1);
	jpgenc_write(state, &len, sizeof(uint16_t), 1);
	jpgenc_write(state, &tc_th, sizeof(uint8_t), 1);
	jpgenc_write(state, matrix_len, sizeof(uint8_t), 16);
	jpgenc_write(state, matrix_val, sizeof(uint8_t), (size_t)num_values);
}


static void jpgenc_encode_and_write_MCU(JPGEncObject *state,
                                        const float *mcu,
                                        JPGEncHuffmanTableType tableType,
                                        int *pred,  // Previous DC coefficient
                                        uint32_t *bitbuffer,  // Bitstack.
                                        uint32_t *location)
{

	uint8_t *huff_dc_len = NULL;
	uint16_t *huff_dc_code = NULL;
	uint8_t *huff_ac_len = NULL;
	uint16_t *huff_ac_code = NULL;

	if (tableType == JPGENC_LUMA) {
		huff_dc_len = state->ehuffsize[JPGENC_LUMA_DC];
		huff_dc_code = state->ehuffcode[JPGENC_LUMA_DC];
		huff_ac_len = state->ehuffsize[JPGENC_LUMA_AC];
		huff_ac_code = state->ehuffcode[JPGENC_LUMA_AC];
	} else {
		huff_dc_len = state->ehuffsize[JPGENC_CHROMA_DC];
		huff_dc_code = state->ehuffcode[JPGENC_CHROMA_DC];
		huff_ac_len = state->ehuffsize[JPGENC_CHROMA_AC];
		huff_ac_code = state->ehuffcode[JPGENC_CHROMA_AC];
	}


	int du[64];  // Data unit in zig-zag order

	float dct_mcu[64];
	memcpy(dct_mcu, mcu, 64 * sizeof(float));

	jpgenc_fdct(dct_mcu);

	for ( int i = 0; i < 64; ++i ) {
		float fval = dct_mcu[i];
		fval *= pqt_matrix[i];
		du[zig_zag_order[i]] = (int) lroundf(fval);
	}

	uint16_t vli[2];

	// Encode DC coefficient.
	int diff = du[0] - *pred;
	*pred = du[0];

	if ( diff != 0 ) {
		jpgenc_calculate_variable_length_int(diff, vli);
		// Write number of bits with Huffman coding
		jpgenc_write_bits(state, bitbuffer, location, huff_dc_len[vli[1]], huff_dc_code[vli[1]]);
		// Write the bits.
		jpgenc_write_bits(state, bitbuffer, location, vli[1], vli[0]);
	} else {
		jpgenc_write_bits(state, bitbuffer, location, huff_dc_len[0], huff_dc_code[0]);
	}

	// ==== Encode AC coefficients ====

	int last_non_zero_i = 0;

	// Find the last non-zero element.
	for ( int i = 63; i > 0; --i ) {
		if (du[i] != 0) {
			last_non_zero_i = i;
			break;
		}
	}

	for ( int i = 1; i <= last_non_zero_i; ++i ) {
		// If zero, increase count. If >=15, encode (FF,00)
		int zero_count = 0;

		while ( du[i] == 0 ) {
			++zero_count;
			++i;

			if (zero_count == 16) {
				// encode (ff,00) == 0xf0
				jpgenc_write_bits(state, bitbuffer, location, huff_ac_len[0xf0], huff_ac_code[0xf0]);
				zero_count = 0;
			}
		}

		jpgenc_calculate_variable_length_int(du[i], vli);

		assert(zero_count < 0x10);
		assert(vli[1] <= 10);

		uint16_t sym1 = (uint16_t)((uint16_t)zero_count << 4) | vli[1];

		assert(huff_ac_len[sym1] != 0);

		// Write symbol 1  --- (RUNLENGTH, SIZE)
		jpgenc_write_bits(state, bitbuffer, location, huff_ac_len[sym1], huff_ac_code[sym1]);
		// Write symbol 2  --- (AMPLITUDE)
		jpgenc_write_bits(state, bitbuffer, location, vli[1], vli[0]);
	}

	if (last_non_zero_i != 63) {
		// write EOB HUFF(00,00)
		jpgenc_write_bits(state, bitbuffer, location, huff_ac_len[0], huff_ac_code[0]);
	}

	return;
}

static JPGEncJPEGHeader create_jpeg_header(void)
{
	JPGEncJPEGHeader header;
	// JFIF header.
	header.SOI = htons(START_MARKER);  // Sequential DCT
	header.APP0 = htons(APP0_MARKER);

	uint16_t jfif_len = sizeof(JPGEncJPEGHeader) - 4 /*SOI & APP0 markers*/;
	header.jfif_len = htons(jfif_len);
	strcpy(header.jfif_id, "JFIF");
	header.version = htons(0x0102);
	header.units = 1;  // Dots-per-inch
	header.x_density = htons(96);
	header.y_density = htons(96);
	header.x_thumb = 0;
	header.y_thumb = 0;
	return header;
}

static JPGEncJPEGComment create_comment_header(void)
{
	JPGEncJPEGComment com;
	uint16_t com_len = 2 + sizeof(jpgenc_com_str) - 1;
	// Comment
	com.com = htons(COMMENT_MARKER);
	com.com_len = htons(com_len);
	memcpy(com.com_str, (const void *) jpgenc_com_str, sizeof(jpgenc_com_str) - 1);
	return com;
}

static JPGEncFrameHeader create_frame_header(int width, int height)
{
	JPGEncFrameHeader header;
	header.SOF = htons(SOF0_MARKER);
	header.len = htons(8 + 3 * 3);
	header.precision = 8;
	assert(width <= UINT16_MAX);
	assert(height <= UINT16_MAX);
	header.width = htons((uint16_t) width);
	header.height = htons((uint16_t) height);
	header.num_components = 3;
	uint8_t tables[3] = {
		0,  // Luma component gets luma table (see jpgenc_write_DQT call above.)
		1,  // Chroma component gets chroma table
		1,  // Chroma component gets chroma table
	};

	for (int i = 0; i < 3; ++i) {
		JPGEncComponentSpec spec;
		spec.component_id = (uint8_t)(i + 1);  // No particular reason. Just 1, 2, 3.
		spec.sampling_factors = (uint8_t) 0x11;
		spec.qt = tables[i];

		header.component_spec[i] = spec;
	}

	return header;
}

static JPGEncScanHeader create_scan_header(void)
{
	JPGEncScanHeader header;
	header.SOS = htons(SOS_MARKER);
	header.len = htons((uint16_t)(6 + (sizeof(JPGEncFrameComponentSpec) * 3)));
	header.num_components = 3;

	uint8_t tables[3] = {
		0x00,
		0x11,
		0x11,
	};

	for (int i = 0; i < 3; ++i) {
		JPGEncFrameComponentSpec cs;
		// Must be equal to component_id from frame header above.
		cs.component_id = (uint8_t)(i + 1);
		cs.dc_ac = (uint8_t)tables[i];

		header.component_spec[i] = cs;
	}

	header.first = 0;
	header.last  = 63;
	header.ah_al = 0;
	return header;
}

static void pixel_encode_rgb_rgba(int src_index, const uint8_t *src_data, float *luma, float *cb, float *cr)
{
	uint8_t r = src_data[src_index + 0];
	uint8_t g = src_data[src_index + 1];
	uint8_t b = src_data[src_index + 2];

	*luma = 0.299f   * r + 0.587f    * g + 0.114f    * b - 128;
	*cb   = -0.1687f * r - 0.3313f   * g + 0.5f      * b;
	*cr   = 0.5f     * r - 0.4187f   * g - 0.0813f   * b;
}

static void pixel_encode_yuyv(int src_index, const uint8_t *src_data, float *luma, float *cb, float *cr)
{
	*luma = (src_data[src_index + 0] - 127);
	int noparity = src_index % 4;

	if (noparity == 0) {
		*cb   = src_data[src_index + 1] - 127;
		*cr   = src_data[src_index + 3] - 127;
	} else {
		*cb   = src_data[src_index - 1] - 127;
		*cr   = src_data[src_index + 1] - 127;
	}
}

static void write_block_encode(JPGEncObject *state, int x, int y, int width, int height,
                               JPGEncFormat src_format, const uint8_t *src_data,
                               int *pred_y, int *pred_b, int *pred_r,
                               uint32_t *bitbuffer, uint32_t *location)
{
	// Block loop: ====
	float du_y[64];
	float du_b[64];
	float du_r[64];

	for ( int off_y = 0; off_y < 8; ++off_y ) {
		for ( int off_x = 0; off_x < 8; ++off_x ) {
			int block_index = (off_y * 8 + off_x);

			float luma = 0, cb = 0, cr = 0;
			int lp = (int) src_format; //количество байн на пиксель

			int src_index = (((y + off_y) * width) + (x + off_x)) * lp;

			int col = x + off_x;
			int row = y + off_y;

			if (row >= height) {
				src_index -= (width * (row - height + 1)) * lp;
			}

			if (col >= width) {
				src_index -= (col - width + 1) * lp;
			}

			assert(src_index < width * height * lp);

			switch (src_format) {
				case JPGEncFormat_RGB:
				case JPGEncFormat_RGBA:
					pixel_encode_rgb_rgba(src_index, src_data, &luma, &cb, &cr);
					break;

				case JPGEncFormat_YUYV:
					pixel_encode_yuyv(src_index, src_data, &luma, &cb, &cr);
					break;
			}

			du_y[block_index] = luma;
			du_b[block_index] = cb;
			du_r[block_index] = cr;
		}
	}

	jpgenc_encode_and_write_MCU(state, du_y, JPGENC_LUMA, pred_y, bitbuffer, location);
	jpgenc_encode_and_write_MCU(state, du_b, JPGENC_CHROMA, pred_b, bitbuffer, location);
	jpgenc_encode_and_write_MCU(state, du_r, JPGENC_CHROMA, pred_r, bitbuffer, location);
}


static int jpgenc_encode_main(JPGEncObject *state,
                              const unsigned char *src_data,
                              int width,
                              int height,
                              JPGEncFormat src_format)
{
	if (width > UINT16_MAX || height > UINT16_MAX) {
		return 0;
	}

	// Write header
	JPGEncJPEGHeader jpeg_header = create_jpeg_header();
	jpgenc_write(state, &jpeg_header, sizeof(JPGEncJPEGHeader), 1);

	// Write comment
	JPGEncJPEGComment com = create_comment_header();
	jpgenc_write(state, &com, sizeof(JPGEncJPEGComment), 1);

	// Write quantization tables.
	jpgenc_write_DQT(state, 0x00);
	jpgenc_write_DQT(state, 0x01);

	// Write the frame marker.
	JPGEncFrameHeader frame_header = create_frame_header(width, height);
	jpgenc_write(state, &frame_header, sizeof(JPGEncFrameHeader), 1);

	jpgenc_write_DHT(state, state->ht_bits[JPGENC_LUMA_DC],   state->ht_vals[JPGENC_LUMA_DC], JPGENC_DC, 0);
	jpgenc_write_DHT(state, state->ht_bits[JPGENC_LUMA_AC],   state->ht_vals[JPGENC_LUMA_AC], JPGENC_AC, 0);
	jpgenc_write_DHT(state, state->ht_bits[JPGENC_CHROMA_DC], state->ht_vals[JPGENC_CHROMA_DC], JPGENC_DC, 1);
	jpgenc_write_DHT(state, state->ht_bits[JPGENC_CHROMA_AC], state->ht_vals[JPGENC_CHROMA_AC], JPGENC_AC, 1);

	// Write start of scan
	JPGEncScanHeader scan_header = create_scan_header();
	jpgenc_write(state, &scan_header, sizeof(JPGEncScanHeader), 1);

	// Set diff to 0.
	int pred_y = 0;
	int pred_b = 0;
	int pred_r = 0;

	// Bit stack
	uint32_t bitbuffer = 0;
	uint32_t location = 0;

	for ( int y = 0; y < height; y += 8 ) {
		for ( int x = 0; x < width; x += 8 ) {
			write_block_encode(state, x, y, width, height, src_format, src_data, &pred_y, &pred_b, &pred_r, &bitbuffer, &location);
		}
	}

	// Finish the image.
	{
		// Flush
		if (location > 0 && location < 8) {
			jpgenc_write_bits(state, &bitbuffer, &location, (uint16_t)(8 - location), 0);
		}
	}
	uint16_t EOI = htons(END_OF_IMAGE_MARKER);
	jpgenc_write(state, &EOI, sizeof(uint16_t), 1);

	if (state->output_buffer_count) {
		state->write_context.func(state->write_context.context, state->output_buffer, state->output_buffer_count);
		state->output_buffer_count = 0;
	}

	return 1;
}

static bool tje_prepare_state(JPGEncObject *state, jpgenc_write_func *func, void *context)
{
	memset(state, 0, sizeof (JPGEncObject));

	JPGEncWriteContext wc = { 0 };

	wc.context = context;
	wc.func = func;

	state->write_context = wc;

	jpgenc_huff_expand(state);
	return true;
}


int jpgenc_encode(jpgenc_write_func *func, void *context, int width, int height, JPGEncFormat src_format, const uint8_t *src_data)
{
	JPGEncObject state;
	tje_prepare_state(&state, func, context);

	int result = jpgenc_encode_main(&state, src_data, width, height, src_format);

	return result;
}
