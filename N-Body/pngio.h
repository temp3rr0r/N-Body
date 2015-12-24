// Copyright (c) 2015, Przemyslaw Klosiewicz
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#ifndef PNGIO_H
#define PNGIO_H

/**
 * This code provides basic reading and writing of PNG files to an easy
 * structure called ImageRGB. The process should be foolproof provided no
 * exotic PNGs are used.
 *
 * libpng is used for IO and the code uses C++11 for some stuff so take the
 * necessary precautions to compile and link it.
 */

#include <cstdint>
#include <string>
#include <vector>
#include <array>

/**
 * This structure defines a simple image encoded as a vector of RGB triples of
 * 8 bits per color. Width is the number of horizontal pixels, while height is
 * the number of vertical pixels. Image data is stored row-wise.
 *
 * Accessing the rgb values for a pixel in (row,column) (i,j) goes as:
 *
 * ImageRGB img = readpng("image.png");
 *
 * uint8_t r = get<0>(img.data[i * img.width + j]);
 * uint8_t g = get<1>(img.data[i * img.width + j]);
 * uint8_t b = get<2>(img.data[i * img.width + j]);
 */
struct ImageRGB {
    uint32_t width;
    uint32_t height;
    std::vector<std::array<uint8_t,3>> data;
};

/**
 * Reads a PNG image from file named 'file_name' and returns the corresponding
 * pixel data as an ImageRGB object.
 * Errors are catched internally but error messages do echo on screen.
 * In case of trouble an 'empty' ImageRGB is returned.
 */
ImageRGB readpng(const std::string & file_name);

/**
 * Saves an ImageRGB object to a PNG file named 'file_name'.
 * Errors are catched internally but error messages do echo on screen.
 */
void writepng(const std::string & file_name, const ImageRGB & image);

#endif /* PNGIO_H */
