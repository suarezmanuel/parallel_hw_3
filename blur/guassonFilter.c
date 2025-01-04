//
//  guassonFilter.c
//  guasonFilter
//
//  Created by Saar Azari on 04/07/2024.
//

#include "guassonFilter.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// Function to load an image
Image *loadImage(const char *filename) {
    int width, height, channels;
    unsigned char *data = stbi_load(filename, &width, &height, &channels, 4); // 4 means we want the image in RGBA format
    if (data == NULL) {
        printf("Failed to load image: %s\n", filename);
        return NULL;
    }

    Image *image = (Image *)malloc(sizeof(Image));
    image->width = width;
    image->height = height;
    image->pixels = (RGBAf *)malloc(width * height * sizeof(RGBAf));

    for (int i = 0; i < width * height; i++) {
        image->pixels[i].r = (float) data[4 * i + 0];
        image->pixels[i].g = (float) data[4 * i + 1];
        image->pixels[i].b = (float) data[4 * i + 2];
        image->pixels[i].a = (float) data[4 * i + 3];
    }

    stbi_image_free(data);
    return image;
}

// Function to save an image
void saveImage(const char *filename, Image *image) {
    unsigned char *data = (unsigned char *)malloc(image->width * image->height * 4);

    for (int i = 0; i < image->width * image->height; i++) {
        data[4 * i + 0] = image->pixels[i].r;
        data[4 * i + 1] = image->pixels[i].g;
        data[4 * i + 2] = image->pixels[i].b;
        data[4 * i + 3] = image->pixels[i].a;
    }

    if (stbi_write_png(filename, image->width, image->height, 4, data, image->width * 4)) {
        printf("Image saved successfully: %s\n", filename);
    } else {
        printf("Failed to save image: %s\n", filename);
    }

    free(data);
}

void createGaussianKernel2D(int radius, double sigma, double **kernel) {
    int kernelWidth = (2 * radius) + 1;
    double sum = 0.0;

    *kernel = (double *)malloc(kernelWidth * kernelWidth * sizeof(double));

    for (int x = -radius; x <= radius; x++) {
        for (int y = -radius; y <= radius; y++) {
            double exponentNumerator = -(x * x + y * y);
            double exponentDenominator = (2 * sigma * sigma);
            double eExpression = exp(exponentNumerator / exponentDenominator);
            double kernelValue = (eExpression / (2 * M_PI * sigma * sigma));
            (*kernel)[(x + radius) * kernelWidth + (y + radius)] = kernelValue;
            sum += kernelValue;
        }
    }

    for (int i = 0; i < kernelWidth * kernelWidth; i++) {
        (*kernel)[i] /= sum;
    }
}


void createGaussianKernel1D(int radius, double sigma, double **kernel) {
    int kernelWidth = (2 * radius) + 1;
    double sum = 0.0;

    *kernel = (double *)malloc(kernelWidth * sizeof(double));

    for (int x = -radius; x <= radius; x++) {
        double exponentNumerator = -(x * x);
        double exponentDenominator = (2 * sigma * sigma);
        double eExpression = exp(exponentNumerator / exponentDenominator);
        double kernelValue = (eExpression / (2 * M_PI * sigma * sigma));
        (*kernel)[x + radius] = kernelValue;
        sum += kernelValue;
    }

    for (int i = 0; i < kernelWidth; i++) {
        (*kernel)[i] /= sum;
    }
}

void createFloatGaussianKernel1D(int radius, float sigma, float **kernel) {
    int kernelWidth = (2 * radius) + 1;
    double sum = 0.0;

    *kernel = (float *)malloc(kernelWidth * sizeof(float));

    for (int x = -radius; x <= radius; x++) {
        double exponentNumerator = -(x * x);
        double exponentDenominator = (2 * sigma * sigma);
        double eExpression = exp(exponentNumerator / exponentDenominator);
        double kernelValue = (eExpression / (2 * M_PI * sigma * sigma));
        (*kernel)[x + radius] = kernelValue;
        sum += kernelValue;
    }

    for (int i = 0; i < kernelWidth; i++) {
        (*kernel)[i] /= sum;
    }
}


void GaussianKernelFloat1DToSIMD(float* in, __m128** out, int w) {
    
    *out = (__m128*) malloc (w * sizeof(__m128));

    for (int i=0; i < w; i++) {
        float v = in[i];
        (*out)[i] = _mm_set_ps(v, v, v, v);
    }
}