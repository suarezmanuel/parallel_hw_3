#include "guassonFilter.h"
#include "blur.h"
#include <stdlib.h>
#include <stdio.h>

int main() {

    // printf("%d", (128/8) / sizeof(float));

    const char *inputFilename = "flowers.png";
    const char *outputFilename1 = "output_image_serial.png";
    const char *outputFilename2 = "output_image_parallel.png";

    int blurRadius = 5;

    Imagef *inputImage = loadImage(inputFilename);
    if (inputImage == NULL) {
        printf("Failed to load image\n");
        return -1;
    }

    printf("test was %s", test(inputImage) ? "sucessful\n" : "unsuccessful\n");

    benchmark(blurRadius, inputImage, serial, "serial");
    benchmark(blurRadius, inputImage, serial_two_passes, "serial_two_passes");
    benchmark(blurRadius, inputImage, serial_t_p_transpose, "serial_t_p_transpose");
    benchmark(blurRadius, inputImage, serial_t_p_unrolling, "serial_t_p_unrolling");
    benchmark(blurRadius, inputImage, serial_t_p_t_unrolling, "serial_t_p_t_unrolling");
    benchmark(blurRadius, inputImage, serial_t_p_t_simd, "serial_t_p_t_simd");
    benchmark(blurRadius, inputImage, serial_t_p_t_u_omp, "serial_t_p_t_u_omp");
    benchmark(blurRadius, inputImage, serial_t_p_t_s_omp, "serial_t_p_t_s_omp");

    // Image *outputImage1 = createBlurredImage(blurRadius, inputImage);

    Image *outputImage2 = createBlurredImageAbs(blurRadius, inputImage, serial_t_p_t_u_omp);

    // printf("the images by serial, parallel are %s\n\n", (equalImages(outputImage1, outputImage2) ? "equal!!" : "different :c "));

    // saveImage(outputFilename1, outputImage1);

    saveImage(outputFilename2, outputImage2);

    free(inputImage->pixels);
    free(inputImage);
    // free(outputImage1->pixels);
    // free(outputImage1);
    free(outputImage2->pixels);
    free(outputImage2);

    return 0;
}
