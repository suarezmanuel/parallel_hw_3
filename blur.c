#include "blur.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void serial(int radius, Image *image, RGBAd* td, RGBAf* tf, int width, int height, double* kernel, int kernelWidth, Image* outputImage) {
    
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
                    RGBA pixel = image->pixels[imageY * width + imageX];

                    redValue += pixel.r * kernelValue;
                    greenValue += pixel.g * kernelValue;
                    blueValue += pixel.b * kernelValue;
                }
            }

            RGBA *outputPixel = &outputImage->pixels[y * width + x];
            outputPixel->r = (unsigned char)redValue;
            outputPixel->g = (unsigned char)greenValue;
            outputPixel->b = (unsigned char)blueValue;
            outputPixel->a = image->pixels[y * width + x].a; // Preserve the alpha channel
        }
    }
}

void serial_two_passes(int r, Image *i, RGBAd* td, RGBAf* tf, int w, int h, double* k, int k_w, Image* o) {

    for (int y = 0; y < h; y++) {

        for (int x = r; x < w - r; x++) {

            double red_val = 0.0, green_val = 0.0, blue_val = 0.0;

            // # pragma omp simd
            for (int k_x = -r; k_x <= r; k_x++) {
                // where we sample from
                int i_x = x + k_x;
                int i_y = y;

                double k_v = k[k_x + r]; 

                // sample
                RGBA p = i->pixels[i_y * w + i_x];

                // convolution
                red_val   += p.r * k_v;
                green_val += p.g * k_v;
                blue_val  += p.b * k_v;
            }

            // output pixel
            RGBAd* o_p = &td[y*w + x];
            // preserve alpha channel
            o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i->pixels[y*w + x].a;
        }
    }

    for (int y = r; y < h - r; y++) {

        for (int x = r; x < w - r; x++) {

            double red_val = 0.0, green_val = 0.0, blue_val = 0.0;

            for (int k_y = -r; k_y <= r; k_y++) {
                // where we sample from
                int i_x = x;
                int i_y = y + k_y;

                // curr val of kernel
                double k_v = k[k_y + r]; 

                // sample
                RGBAd p = td[i_y * w + i_x];

                // convolution
                red_val   += p.r * k_v;
                green_val += p.g * k_v;
                blue_val  += p.b * k_v;
            }

            // output pixel
            RGBA* o_p = &o->pixels[y*w + x];
            // preserve alpha channel
            o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i->pixels[y*w + x].a;
        }
    }
}

void serial_t_p_transpose(int r, Image *i, RGBAd* td, RGBAf* tf, int w, int h, double* k, int k_w, Image* o) {

    for (int y = 0; y < h; y++) {

        for (int x = r; x < w - r; x++) {

            double red_val = 0.0, green_val = 0.0, blue_val = 0.0;

            // # pragma omp simd
            for (int k_x = -r; k_x <= r; k_x++) {
                // where we sample from
                int i_x = x + k_x;
                int i_y = y;

                double k_v = k[k_x + r]; 

                // sample
                RGBA p = i->pixels[i_y * w + i_x];

                // convolution
                red_val   += p.r * k_v;
                green_val += p.g * k_v;
                blue_val  += p.b * k_v;
            }

            // output pixel
            RGBAd* o_p = &td[x*h + y];
            // preserve alpha channel
            o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i->pixels[y*w + x].a;
        }
    }

    for (int y = r; y < h - r; y++) {

        for (int x = r; x < w - r; x++) {

            double red_val = 0.0, green_val = 0.0, blue_val = 0.0;

            for (int k_y = -r; k_y <= r; k_y++) {
                // where we sample from
                int i_x = x;
                int i_y = y + k_y;

                // curr val of kernel
                double k_v = k[k_y + r]; 

                // sample
                RGBAd p = td[i_x * h + i_y];

                // convolution
                red_val   += p.r * k_v;
                green_val += p.g * k_v;
                blue_val  += p.b * k_v;
            }

            // output pixel
            RGBA* o_p = &o->pixels[y*w + x];
            // preserve alpha channel
            o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i->pixels[y*w + x].a;
        }
    }
}

void serial_t_p_unrolling(int r, Image *i, RGBAd* td, RGBAf* tf, int w, int h, double* k, int k_w, Image* o) {

    for (int y = 0; y < h; y++) {
        for (int x = r; x < w - r; x++) {

            float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

            int kx;

            for (kx = -r; kx <= r-3; kx += 4) {
                float k_v0 = k[kx     + r];
                float k_v1 = k[kx + 1 + r];
                float k_v2 = k[kx + 2 + r];
                float k_v3 = k[kx + 3 + r];

                RGBA p0 = i->pixels[y * w + (x + kx)];
                RGBA p1 = i->pixels[y * w + (x + kx + 1)];
                RGBA p2 = i->pixels[y * w + (x + kx + 2)];
                RGBA p3 = i->pixels[y * w + (x + kx + 3)];

                red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
                green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
                blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
            }

            for (; kx <= r; kx++) {
                float k_v = k[kx + r];
                RGBA p = i->pixels[y * w + (x + kx)];
                red_val   += p.r * k_v;
                green_val += p.g * k_v;
                blue_val  += p.b * k_v;
            }

            RGBAd *o_p = &td[y * w + x]; o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i->pixels[y * w + x].a;
        }
    }

    for (int y = r; y < h - r; y++) {
        for (int x = r; x < w - r; x++) {

            float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

            int ky;

            for (ky = -r; ky <= r-3; ky += 4) {

                float k_v0 = k[ky + r];
                float k_v1 = k[ky + r + 1];
                float k_v2 = k[ky + r + 2];
                float k_v3 = k[ky + r + 3];

                RGBAd p0 = td[(y + ky) * w + x];
                RGBAd p1 = td[(y + ky + 1) * w + x];
                RGBAd p2 = td[(y + ky + 2) * w + x];
                RGBAd p3 = td[(y + ky + 3) * w + x];

                red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
                green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
                blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
            }

            for (; ky <= r; ky++) {
                float k_v = k[ky + r];
                RGBAd p = td[(y + ky) * w + x]; 
                red_val   += p.r * k_v;
                green_val += p.g * k_v;
                blue_val  += p.b * k_v;
            }

            RGBA* out = &o->pixels[y * w + x]; out->r = (unsigned char)red_val; out->g = (unsigned char)green_val; out->b = (unsigned char)blue_val; out->a = i->pixels[y * w + x].a;
        }
    }
}

void serial_t_p_t_simd(int r, Image *i, RGBAd* td, RGBAf* tf, int w, int h, double* k, int k_w, Image* o) {

    for (int y = 0; y < h; y++) {
        for (int x = r; x < w - r; x++) {

            float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

            int kx;

            for (kx = -r; kx <= r-3; kx += 4) {
                float k_v0 = k[kx     + r];
                float k_v1 = k[kx + 1 + r];
                float k_v2 = k[kx + 2 + r];
                float k_v3 = k[kx + 3 + r];

                RGBA p0 = i->pixels[y * w + (x + kx)];
                RGBA p1 = i->pixels[y * w + (x + kx + 1)];
                RGBA p2 = i->pixels[y * w + (x + kx + 2)];
                RGBA p3 = i->pixels[y * w + (x + kx + 3)];

                red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
                green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
                blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
            }

            for (; kx <= r; kx++) {
                float k_v = k[kx + r];
                RGBA p = i->pixels[y * w + (x + kx)];
                red_val   += p.r * k_v;
                green_val += p.g * k_v;
                blue_val  += p.b * k_v;
            }

            RGBAf *o_p = &tf[x * h + y]; o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i->pixels[y * w + x].a;
        }
    }

    for (int y = r; y < h - r; y++) {
        for (int x = r; x < w - r; x++) {

            float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

            int ky;

            for (ky = -r; ky <= r-3; ky += 4) {

                float k_v0 = k[ky + r];
                float k_v1 = k[ky + r + 1];
                float k_v2 = k[ky + r + 2];
                float k_v3 = k[ky + r + 3];

                RGBAf p0 = tf[x * h + y + ky];
                RGBAf p1 = tf[x * h + y + ky + 1];
                RGBAf p2 = tf[x * h + y + ky + 2];
                RGBAf p3 = tf[x * h + y + ky + 3];

                red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
                green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
                blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
            }

            for (; ky <= r; ky++) {
                float k_v = k[ky + r];
                RGBAf p = tf[x * h + y + ky]; 
                red_val   += p.r * k_v;
                green_val += p.g * k_v;
                blue_val  += p.b * k_v;
            }

            RGBA* out = &o->pixels[y * w + x]; out->r = (unsigned char)red_val; out->g = (unsigned char)green_val; out->b = (unsigned char)blue_val; out->a = i->pixels[y * w + x].a;
        }
    }
}


Image *createBlurredImageAbs(int radius, Image *image, blur_func f) {

    int width = image->width;
    int height = image->height;
    Image *outputImage = (Image *)malloc(sizeof(Image));
    outputImage->width = width;
    outputImage->height = height;
    outputImage->pixels = (RGBA *)calloc(width * height, sizeof(RGBA));

    // standard deviation
    double sigma = fmax(radius / 2.0, 1.0);
    int kernelWidth = (2 * radius) + 1;
    double *kernel1;
    double *kernel2;
    double sum;
    
    RGBAd* td = (RGBAd*) malloc(width * height * sizeof(RGBAd));
    RGBAf* tf = (RGBAf*) malloc(width * height * sizeof(RGBAf));
    
    createGaussianKernel1D(radius, sigma, &kernel1, &sum);
    createGaussianKernel2D(radius, sigma, &kernel2, &sum);

    double* kernel = (f == serial) ? kernel2 : kernel1;
    
    f(radius, image, td, tf, width, height, kernel, kernelWidth, outputImage);

    free(kernel);

    free(td);
    free(tf);

    return outputImage;
}

Image *createBlurredImage(int radius, Image *image) { return createBlurredImageAbs(radius, image, serial); }

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

void benchmark (int radius, Image *image, blur_func f, char* name) {

    double s, e, t;

    for (int i=0; i < BENCHMARK_ITERS; i++) {
        s = (double) clock() / CLOCKS_PER_SEC;

        createBlurredImageAbs(radius, image, f);

        e = (double) clock() / CLOCKS_PER_SEC; 

        t += (e - s) / BENCHMARK_ITERS;
    }

    printf("%-30s took %f seconds\n", name, e-s); 
}