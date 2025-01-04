#include "guassonFilter.h"

#define BENCHMARK_ITERS 100

typedef void (*blur_func) (int, RGBAf*, RGBAd*, RGBAf*, int, int, double*, float*, __m128*, int, RGBAf* o);

Image *createBlurredImage(int radius, Image *image);

Image *createBlurredImageAbs(int radius, Image *image, blur_func f);

void serial(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBAf* o);

void parallel(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBAf* o);

void benchmark (int radius, Image *image, blur_func f, char* name);

int equalImages (Image* i1, Image* i2);