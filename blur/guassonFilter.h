// 340858935 Manuel Ignacio Suarez Del Solar
// 216270827 Azmon Avraham Zvi

//
//  guassonFilter.h
//  guasonFilter
//
//  Created by Saar Azari on 04/07/2024.
//

#ifndef guassonFilter_h
#define guassonFilter_h

// #include <smmintrin.h>
#include <emmintrin.h>
// #include <xmmintrin.h>
// #include <intrin.h>

// Define an RGBA pixel structure
typedef struct {
    unsigned char r, g, b, a;
} RGBA;

typedef union {
    __m128i sse; // first 4 bytes are first RGBA
    unsigned char E[16];
} RGBA_x4;

typedef struct {
    double r, g, b, a;
} RGBAd;

typedef struct {
    float r, g, b, a;
} RGBAf;

typedef union {
    __m128 sse;
    float E [4];
} RGBAf_x1;

// Define an Image structure
// typedef struct {
//     int width;
//     int height;
//     RGBA *pixels;
// } Image;

typedef struct {
    int width;
    int height;
    RGBAf *pixels;
} Image;

Image *loadImage(const char *filename);

void saveImage(const char *filename, Image *image);

void createGaussianKernel1D(int radius, double sigma, double **kernel);

void createFloatGaussianKernel1D(int radius, float sigma, float **kernel);

void GaussianKernelFloat1DToSIMD(float* in, __m128** out, int w);

void createGaussianKernel2D(int radius, double sigma, double **kernel);

#endif /* guassonFilter_h */
