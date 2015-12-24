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
#include "pngio.h"
#include <png.h>
#include <iostream>
#include <stdexcept>

using namespace std;

ImageRGB readpng(const string & file_name) {
    // PNG structure & info containers
    png_structp png = NULL;
    png_infop info = NULL;
    // Raw file pointer used by libpng
    FILE * fp = NULL;
    // Data will be stored in an ImageRGB struct object
    ImageRGB image;

    try {
        // Open file for reading
        fp = fopen(file_name.c_str(), "rb");
        if (!fp) throw runtime_error("File \'" + file_name + "\' not found");

        // 8 byte PNG signature check
        unsigned char sig[8];
        size_t _ = fread(sig, 1, 8, fp);
        if (!png_check_sig(sig, 8)) throw runtime_error("Bad PNG signature");

        // Allocate PNG read and info structures with default error handlers
        png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png) { throw runtime_error("Out of memory"); }
        info = png_create_info_struct(png);
        if (!info) { throw runtime_error("Out of memory"); }

        // Old-school error handler with setjmp/longjmp. I don't know if
        // throwing an exception here is 100% safe but I don't really care.
        if (setjmp(png_jmpbuf(png))) { throw runtime_error("Error"); }

        // Initialize default IO functions
        png_init_io(png, fp);

        // Notify libpng we already read 8 bytes from fp to check the PNG sig.
        png_set_sig_bytes(png, 8);

        // Use high-level reading & ensure 24bpp with some basic transforms
        png_read_png(png, info, PNG_TRANSFORM_STRIP_16 |
                                PNG_TRANSFORM_STRIP_ALPHA, NULL);
        png_bytep * row_pointers = png_get_rows(png, info);

        // Get & store image dimensions
        png_uint_32 width = png_get_image_width(png, info);
        png_uint_32 height = png_get_image_height(png, info);
        image.width = width;
        image.height = height;

        // Ready to transform data to RGB triples
        for (size_t i = 0; i < height; ++i) {
            for (size_t j = 0; j < width*3; j += 3) {
                const uint8_t r = row_pointers[i][j];
                const uint8_t g = row_pointers[i][j+1];
                const uint8_t b = row_pointers[i][j+2];
                image.data.push_back({r,g,b});
            }
        }

        // All went good, clean up
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
    } catch (exception & e) {
        cerr << "In function \'" << __func__ << "\': " << e.what() << endl;
        // Attempt some form of graceful cleanup
        if (png) { png_destroy_read_struct(&png, &info, NULL); }
        if (fp) { fclose(fp); }
        throw;
    }
    // Get the hell out, with or without a proper image
    return image;
}

void writepng(const string & file_name, const ImageRGB & image) {
    // PNG structure & info containers
    png_structp png = NULL;
    png_infop info = NULL;
    // Raw file pointer used by libpng
    FILE * fp = NULL;

    try {
        // Basic sanity checks
        if (!image.data.size()) { throw runtime_error("Empty image data"); }
        if (file_name.size() == 0) { throw runtime_error("Empty file_name"); }
        if (image.data.size() != image.width * image.height) {
            throw runtime_error("Img size doesn\'t match number of pixels");
        }

        // Open file for writing
        fp = fopen(file_name.c_str(), "wb");
        if (!fp) throw runtime_error("Can\'t write to \'" + file_name + "\'");

        // Allocate PNG write and info structures with default error handlers
        png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png) { throw runtime_error("Out of memory"); }
        info = png_create_info_struct(png);
        if (!info) { throw runtime_error("Out of memory"); }

        // Old-school error handler with setjmp/longjmp. I don't know if
        // throwing an exception here is 100% safe but I don't really care.
        if (setjmp(png_jmpbuf(png))) { throw runtime_error("Error"); }

        // Initialize default IO functions
        png_init_io(png, fp);

        // Fill the info structure with image metadata
        png_set_IHDR(png, info, image.width, image.height,
                     8, PNG_COLOR_TYPE_RGB, // 8bit per pixel RGB
                     PNG_INTERLACE_NONE,    // No interlacing
                     PNG_COMPRESSION_TYPE_DEFAULT,
                     PNG_FILTER_TYPE_DEFAULT);

        png_bytep * row_pointers = new png_bytep[image.height];
        for (size_t i = 0; i < image.height; ++i) {
            row_pointers[i] = new png_byte[3*image.width];
            for (size_t j = 0; j < image.width; ++j) {
                const size_t k = i * image.width + j;
                row_pointers[i][3*j] = get<0>(image.data[k]);
                row_pointers[i][3*j+1] = get<1>(image.data[k]);
                row_pointers[i][3*j+2] = get<2>(image.data[k]);
            }
        }
        png_set_rows(png, info, row_pointers);

        // Use high-level writing
        png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);

        // All went good, clean up & free memory
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        for (size_t i = 0; i < image.height; ++i) { delete[] row_pointers[i]; }
        delete[] row_pointers;
    } catch (exception & e) {
        cerr << "In function \'" << __func__ << "\': " << e.what() << endl;
        // Attempt some form of graceful cleanup
        // NOTE: I don't free memory allocated for row_pointers: memleak!
        if (png) { png_destroy_write_struct(&png, &info); }
        if (fp) { fclose(fp); }
        throw;
    }
}

