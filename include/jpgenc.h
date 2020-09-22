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

#ifndef JPGENC_H
#define JPGENC_H

#include <stdint.h>
#include <stdio.h>

typedef void jpgenc_write_func(void *context, const void *data, size_t size);

typedef enum {
	JPGEncFormat_RGB = 3,
	JPGEncFormat_RGBA = 4,
	JPGEncFormat_YUYV = 2
} JPGEncFormat;

int jpgenc_encode(jpgenc_write_func *func,
                  void *context,
                  int width,
                  int height,
                  JPGEncFormat src_format,
                  const uint8_t *src_data);

#endif // JPGENC_H
