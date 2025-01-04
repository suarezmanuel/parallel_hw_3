#include "blur.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void serial(int radius, RGBAf* pixels, RGBAd* td, RGBAf* tf, int width, int height, double* kernel, float* kf, __m128* k_vv, int kernelWidth, RGBAf* outputPixels) {
    
    for (int x = radius; x < width - radius; x++) {
        for (int y = radius; y < height - radius; y++) {
            double redValue = 0.0;
            double greenValue = 0.0;
            double blueValue = 0.0;

            for (int kernelX = -radius; kernelX <= radius; kernelX++) {
                for (int kernelY = -radius; kernelY <= radius; kernelY++) {
                    int imageX = x - kernelX;
                    int imageY = y - kernelY;
                    double kernelValue = kernel[(kernelX + radius) * kernelWidth + (kernelY + radius)];
                    RGBAf pixel = pixels[imageY * width + imageX];

                    redValue += pixel.r * kernelValue;
                    greenValue += pixel.g * kernelValue;
                    blueValue += pixel.b * kernelValue;
                }
            }

            RGBAf *outputPixel = &outputPixels[y * width + x];
            outputPixel->r = (unsigned char)redValue;
            outputPixel->g = (unsigned char)greenValue;
            outputPixel->b = (unsigned char)blueValue;
            outputPixel->a = pixels[y * width + x].a; // Preserve the alpha channel
        }
    }
}

void parallel(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBAf* o) {

    __m128 _0 = _mm_set_ps1(0);

    #pragma omp parallel 
    {

        #pragma omp for
        for (int y = 0; y < h; y++) {
            for (int x = r; x < w - r; x++) {

                RGBAf_x1 acc; acc.sse = _0;
                int kx;

                // #pragma omp simd
                for (kx = -r; kx <= r-3; kx += 4) {
                    __m128 k_vv0 = k_vv[kx     + r];
                    __m128 k_vv1 = k_vv[kx + 1 + r];
                    __m128 k_vv2 = k_vv[kx + 2 + r];
                    __m128 k_vv3 = k_vv[kx + 3 + r];

                    __m128 p0 = _mm_loadu_ps((float*) &i[y * w + (x + kx)]);
                    __m128 p1 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 1)]);
                    __m128 p2 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 2)]);
                    __m128 p3 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 3)]);

                    p0 = _mm_mul_ps(p0, k_vv0);
                    p1 = _mm_mul_ps(p1, k_vv1);
                    p2 = _mm_mul_ps(p2, k_vv2);
                    p3 = _mm_mul_ps(p3, k_vv3);

                    __m128 s1 = _mm_add_ps(p0, p1);
                    __m128 s2 = _mm_add_ps(p2, p3);

                    acc.sse = _mm_add_ps(acc.sse, s1);
                    acc.sse = _mm_add_ps(acc.sse, s2);
                }

                for (; kx <= r; kx++) {
                    float k_v = kf[kx + r];
                    RGBAf p = i[y * w + (x + kx)];
                    acc.E[0] += p.r * k_v;
                    acc.E[1] += p.g * k_v;
                    acc.E[2] += p.b * k_v;
                }

                // maybe we could even make a temp of vectors, so we don't have to convert
                RGBAf *o_p = &tf[x * h + y]; o_p->r = acc.E[0]; o_p->g = acc.E[1]; o_p->b = acc.E[2]; o_p->a = i[y * w + x].a;
            }
        }

        #pragma omp for
        for (int y = r; y < h - r; y++) {
            for (int x = r; x < w - r; x++) {

                RGBAf_x1 acc; acc.sse = _0;
                int ky;

                // #pragma omp simd
                for (ky = -r; ky <= r-3; ky += 4) {
                    __m128 k_vv0 = k_vv[ky + r    ];
                    __m128 k_vv1 = k_vv[ky + r + 1];
                    __m128 k_vv2 = k_vv[ky + r + 2];
                    __m128 k_vv3 = k_vv[ky + r + 3];

                    __m128 p0 = _mm_loadu_ps((float*) &tf[x * h + (y + ky)]);
                    __m128 p1 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 1)]);
                    __m128 p2 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 2)]);
                    __m128 p3 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 3)]);

                    p0 = _mm_mul_ps(p0, k_vv0);
                    p1 = _mm_mul_ps(p1, k_vv1);
                    p2 = _mm_mul_ps(p2, k_vv2);
                    p3 = _mm_mul_ps(p3, k_vv3);

                    __m128 s1 = _mm_add_ps(p0, p1);
                    __m128 s2 = _mm_add_ps(p2, p3);

                    acc.sse = _mm_add_ps(acc.sse, s1);
                    acc.sse = _mm_add_ps(acc.sse, s2);
                }

                for (; ky <= r; ky++) {
                    float k_v = kf[ky + r];
                    RGBAf p = tf[x * h + (y + ky)];
                    acc.E[0] += p.r * k_v;
                    acc.E[1] += p.g * k_v;
                    acc.E[2] += p.b * k_v;
                }

                RGBAf* out = &o[y * w + x]; out->r = (unsigned char)acc.E[0]; out->g = (unsigned char)acc.E[1]; out->b = (unsigned char)acc.E[2]; out->a = i[y * w + x].a;
            }
        }
    }
}

Image *createBlurredImageAbs(int radius, Image* image, blur_func f) {

    omp_set_num_threads(32);

    int width = image->width;
    int height = image->height;
    Image *outputImage = (Image *)malloc(sizeof(Image));
    outputImage->width = width;
    outputImage->height = height;
    outputImage->pixels = (RGBAf *)calloc(width * height, sizeof(RGBAf));

    // standard deviation
    double sigma = fmax(radius / 2.0, 1.0);
    int kernelWidth = (2 * radius) + 1;
    double *kernel1;
    __m128 *kernel1_vv;
    double *kernel2;
    float *kf;
    
    RGBAd* td = (RGBAd*) malloc(width * height * sizeof(RGBAd));
    RGBAf* tf = (RGBAf*) malloc(width * height * sizeof(RGBAf));
    
    createGaussianKernel1D(radius, sigma, &kernel1);
    createFloatGaussianKernel1D(radius, sigma, &kf);
    GaussianKernelFloat1DToSIMD(kf, &kernel1_vv, kernelWidth);
    createGaussianKernel2D(radius, sigma, &kernel2);

    double* kernel = (f == serial) ? kernel2 : kernel1;
    
    f(radius, image->pixels, td, tf, width, height, kernel, kf, kernel1_vv, kernelWidth, outputImage->pixels);

    free(kernel);
    free(kf);
    free(td);
    free(tf);

    return outputImage;
}

Image *createBlurredImage(int radius, Image* image) { return createBlurredImageAbs(radius, image, parallel); }

int equalImages (Image* i1, Image* i2) {
    for (int i=0; i < i1->width; i++) {
        for (int j=0; j < i1->height; j++) {
            if ((i1->pixels[i * i1->height + j].r - i2->pixels[i * i1->height + j].r) > 1 ||
                (i1->pixels[i * i1->height + j].g - i2->pixels[i * i1->height + j].g) > 1 ||
                (i1->pixels[i * i1->height + j].b - i2->pixels[i * i1->height + j].b) > 1 ||
                (i1->pixels[i * i1->height + j].a - i2->pixels[i * i1->height + j].a) > 1 ) 
            { 
                printf("different on %d %d, r: %d %d g: %d %d b: %d %d a: %d %d\n", i, j,
                       i1->pixels[i * i1->height + j].r,
                       i2->pixels[i * i1->height + j].r,
                       i1->pixels[i * i1->height + j].g,
                       i2->pixels[i * i1->height + j].g,
                       i1->pixels[i * i1->height + j].b,
                       i2->pixels[i * i1->height + j].b,
                       i1->pixels[i * i1->height + j].a,
                       i2->pixels[i * i1->height + j].a); 
                return 0;
            }
        }
    }
    return 1;
}

void benchmark (int radius, Image* image, blur_func f, char* name) {

    double s, e, t;

    for (int i=0; i < BENCHMARK_ITERS; i++) {
        s = (double) omp_get_wtime();

        createBlurredImageAbs(radius, image, f);

        e = (double) omp_get_wtime(); 

        t += (e - s) / BENCHMARK_ITERS;
    }

    printf("%-30s took %f seconds\n", name, e-s); 
}