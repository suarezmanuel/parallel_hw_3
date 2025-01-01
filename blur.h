#include "guassonFilter.h"

#define BENCHMARK_ITERS 100

typedef void (*blur_func) (int, Image*, RGBAd*, RGBAf*, int, int, double*, int, Image*);

Image *createBlurredImage(int radius, Image *image);

Image *createBlurredImageAbs(int radius, Image *image, blur_func f);

void serial(int radius, Image *image, RGBAd* td, RGBAf* tf, int width, int height, double* kernel, int kernelWidth, Image* outputImage);

void serial_two_passes(int radius, Image *image, RGBAd* td, RGBAf* tf, int width, int height, double* kernel, int kernelWidth, Image* outputImage);

void serial_t_p_transpose(int radius, Image *image, RGBAd* td, RGBAf* tf, int width, int height, double* kernel, int kernelWidth, Image* outputImage);

void serial_t_p_unrolling(int radius, Image *image, RGBAd* td, RGBAf* tf, int width, int height, double* kernel, int kernelWidth, Image* outputImage);

void serial_t_p_t_simd(int radius, Image *image, RGBAd* td, RGBAf* tf, int width, int height, double* kernel, int kernelWidth, Image* outputImage);

void benchmark (int radius, Image *image, blur_func f, char* name);

int equalImages (Image* i1, Image* i2);