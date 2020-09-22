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

#ifndef JPEGENC_DEFINITION_H
#define JPEGENC_DEFINITION_H

#include <stdint.h>
#include <stdlib.h>
#include "jpgenc_const.h"

enum {
	JPGENC_LUMA_DC,
	JPGENC_LUMA_AC,
	JPGENC_CHROMA_DC,
	JPGENC_CHROMA_AC,
};

typedef enum {
	JPGENC_DC = 0,
	JPGENC_AC = 1
} JPGEncHuffmanTableClass;

typedef enum {
	JPGENC_LUMA,
	JPGENC_CHROMA
} JPGEncHuffmanTableType;

typedef struct {
	void *context;
	jpgenc_write_func *func;
} JPGEncWriteContext;

typedef struct {
	// Huffman data.
	uint8_t         ehuffsize[4][257];
	uint16_t        ehuffcode[4][256];
	uint8_t const *ht_bits[4];
	uint8_t const *ht_vals[4];

	// fwrite by default. User-defined when using tje_encode_with_func.
	JPGEncWriteContext write_context;

	// Buffered output. Big performance win when using the usual stdlib implementations.
	size_t          output_buffer_count;
	uint8_t         output_buffer[JPGENC_BUFFER_SIZE];
} JPGEncObject;

#pragma pack(push)
#pragma pack(1)
typedef struct {
	uint16_t SOI;
	uint16_t APP0;
	uint16_t jfif_len;
	char  jfif_id[5];
	uint16_t version;
	uint8_t  units;
	uint16_t x_density;
	uint16_t y_density;
	uint8_t  x_thumb;
	uint8_t  y_thumb;
} JPGEncJPEGHeader;

typedef struct {
	uint16_t com;
	uint16_t com_len;
	char     com_str[sizeof(jpgenc_com_str) - 1];
} JPGEncJPEGComment;

// Helper struct for TJEFrameHeader (below).
typedef struct {
	uint8_t  component_id;
	uint8_t  sampling_factors;    // most significant 4 bits: horizontal. 4 LSB: vertical (A.1.1)
	uint8_t  qt;                  // Quantization table selector.
} JPGEncComponentSpec;

typedef struct {
	uint16_t         SOF;
	uint16_t         len;                   // 8 + 3 * frame.num_components
	uint8_t          precision;             // Sample precision (bits per sample).
	uint16_t         height;
	uint16_t         width;
	uint8_t          num_components;        // For this implementation, will be equal to 3.
	JPGEncComponentSpec component_spec[3];
} JPGEncFrameHeader;

typedef struct {
	uint8_t component_id;                 // Just as with TJEComponentSpec
	uint8_t dc_ac;                        // (dc|ac)
} JPGEncFrameComponentSpec;

typedef struct {
	uint16_t              SOS;
	uint16_t              len;
	uint8_t               num_components;  // 3.
	JPGEncFrameComponentSpec component_spec[3];
	uint8_t               first;  // 0
	uint8_t               last;  // 63
	uint8_t               ah_al;  // o
} JPGEncScanHeader;
#pragma pack(pop)


#endif // JPEGENC_DEFINITION_H
