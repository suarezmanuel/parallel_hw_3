//
//  guassonFilter.h
//  guasonFilter
//
//  Created by Saar Azari on 04/07/2024.
//

#ifndef guassonFilter_h
#define guassonFilter_h

// #include <smmintrin.h>
// #include <emmintrin.h>
// #include <xmmintrin.h>
// #include <intrin.h>

// Define an RGBA pixel structure
typedef struct {
    unsigned char r, g, b, a;
} RGBA;

typedef union {
    __m128i sse;
    unsigned char E[16];
} RGBA_x4;

typedef struct {
    double r, g, b, a;
} RGBAd;

typedef struct {
    float r, g, b, a;
} RGBAf;

// Define an Image structure
typedef struct {
    int width;
    int height;
    RGBA *pixels;
} Image;

Image *loadImage(const char *filename);

void saveImage(const char *filename, Image *image);

void createGaussianKernel1D(int radius, double sigma, double **kernel, double *sum);

void createGaussianKernel2D(int radius, double sigma, double **kernel, double *sum);

#endif /* guassonFilter_h */
