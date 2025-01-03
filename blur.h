#include "guassonFilter.h"

#define BENCHMARK_ITERS 100

typedef void (*blur_func) (int, RGBAf*, RGBAd*, RGBAf*, int, int, double*, float*, __m128*, int, RGBA* o);

Image *createBlurredImage(int radius, Imagef *image);

Image *createBlurredImageAbs(int radius, Imagef *image, blur_func f);

void serial(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o);
            
void serial_two_passes(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o);

void serial_t_p_transpose(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o);

void serial_t_p_unrolling(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o);

void serial_t_p_t_unrolling(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o);

void serial_t_p_t_u_omp(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o);

void serial_t_p_t_simd(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o);

void serial_t_p_t_s_omp(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o);

void benchmark (int radius, Imagef *image, blur_func f, char* name);

int equalImages (Image* i1, Image* i2);

int test (Imagef* image);