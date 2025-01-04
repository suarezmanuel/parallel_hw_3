// #include "guassonFilter.h"

// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <omp.h>

// void serial_two_passes(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o) {

//     for (int y = 0; y < h; y++) {

//         for (int x = r; x < w - r; x++) {

//             double red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//             // # pragma omp simd
//             for (int k_x = -r; k_x <= r; k_x++) {
//                 // where we sample from
//                 int i_x = x + k_x;
//                 int i_y = y;

//                 double k_v = k[k_x + r]; 

//                 // sample
//                 RGBAf p = i[i_y * w + i_x];

//                 // convolution
//                 red_val   += p.r * k_v;
//                 green_val += p.g * k_v;
//                 blue_val  += p.b * k_v;
//             }

//             // output pixel
//             RGBAd* o_p = &td[y*w + x];
//             // preserve alpha channel
//             o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i[y*w + x].a;
//         }
//     }

//     for (int y = r; y < h - r; y++) {

//         for (int x = r; x < w - r; x++) {

//             double red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//             for (int k_y = -r; k_y <= r; k_y++) {
//                 // where we sample from
//                 int i_x = x;
//                 int i_y = y + k_y;

//                 // curr val of kernel
//                 double k_v = k[k_y + r]; 

//                 // sample
//                 RGBAd p = td[i_y * w + i_x];

//                 // convolution
//                 red_val   += p.r * k_v;
//                 green_val += p.g * k_v;
//                 blue_val  += p.b * k_v;
//             }

//             // output pixel
//             RGBA* o_p = &o[y*w + x];
//             // preserve alpha channel
//             o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i[y*w + x].a;
//         }
//     }
// }

// void serial_t_p_transpose(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o) {

//     for (int y = 0; y < h; y++) {

//         for (int x = r; x < w - r; x++) {

//             double red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//             // # pragma omp simd
//             for (int k_x = -r; k_x <= r; k_x++) {
//                 // where we sample from
//                 int i_x = x + k_x;
//                 int i_y = y;

//                 double k_v = k[k_x + r]; 

//                 // sample
//                 RGBAf p = i[i_y * w + i_x];

//                 // convolution
//                 red_val   += p.r * k_v;
//                 green_val += p.g * k_v;
//                 blue_val  += p.b * k_v;
//             }

//             // output pixel
//             RGBAd* o_p = &td[x*h + y];
//             // preserve alpha channel
//             o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i[y*w + x].a;
//         }
//     }

//     for (int y = r; y < h - r; y++) {

//         for (int x = r; x < w - r; x++) {

//             double red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//             for (int k_y = -r; k_y <= r; k_y++) {
//                 // where we sample from
//                 int i_x = x;
//                 int i_y = y + k_y;

//                 // curr val of kernel
//                 double k_v = k[k_y + r]; 

//                 // sample
//                 RGBAd p = td[i_x * h + i_y];

//                 // convolution
//                 red_val   += p.r * k_v;
//                 green_val += p.g * k_v;
//                 blue_val  += p.b * k_v;
//             }

//             // output pixel
//             RGBA* o_p = &o[y*w + x];
//             // preserve alpha channel
//             o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i[y*w + x].a;
//         }
//     }
// }

// void serial_t_p_unrolling(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o) {

//     for (int y = 0; y < h; y++) {
//         for (int x = r; x < w - r; x++) {

//             float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//             int kx;

//             for (kx = -r; kx <= r-3; kx += 4) {
//                 float k_v0 = k[kx     + r];
//                 float k_v1 = k[kx + 1 + r];
//                 float k_v2 = k[kx + 2 + r];
//                 float k_v3 = k[kx + 3 + r];

//                 RGBAf p0 = i[y * w + (x + kx)];
//                 RGBAf p1 = i[y * w + (x + kx + 1)];
//                 RGBAf p2 = i[y * w + (x + kx + 2)];
//                 RGBAf p3 = i[y * w + (x + kx + 3)];

//                 red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
//                 green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
//                 blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
//             }

//             for (; kx <= r; kx++) {
//                 float k_v = k[kx + r];
//                 RGBAf p = i[y * w + (x + kx)];
//                 red_val   += p.r * k_v;
//                 green_val += p.g * k_v;
//                 blue_val  += p.b * k_v;
//             }

//             RGBAd *o_p = &td[y * w + x]; o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i[y * w + x].a;
//         }
//     }

//     for (int y = r; y < h - r; y++) {
//         for (int x = r; x < w - r; x++) {

//             float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//             int ky;

//             for (ky = -r; ky <= r-3; ky += 4) {

//                 float k_v0 = k[ky + r];
//                 float k_v1 = k[ky + r + 1];
//                 float k_v2 = k[ky + r + 2];
//                 float k_v3 = k[ky + r + 3];

//                 RGBAd p0 = td[(y + ky) * w + x];
//                 RGBAd p1 = td[(y + ky + 1) * w + x];
//                 RGBAd p2 = td[(y + ky + 2) * w + x];
//                 RGBAd p3 = td[(y + ky + 3) * w + x];

//                 red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
//                 green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
//                 blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
//             }

//             for (; ky <= r; ky++) {
//                 float k_v = k[ky + r];
//                 RGBAd p = td[(y + ky) * w + x]; 
//                 red_val   += p.r * k_v;
//                 green_val += p.g * k_v;
//                 blue_val  += p.b * k_v;
//             }

//             RGBA* out = &o[y * w + x]; out->r = (unsigned char)red_val; out->g = (unsigned char)green_val; out->b = (unsigned char)blue_val; out->a = i[y * w + x].a;
//         }
//     }
// }

// void serial_t_p_t_unrolling(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o) {

//     for (int y = 0; y < h; y++) {
//         for (int x = r; x < w - r; x++) {

//             float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//             int kx;

//             for (kx = -r; kx <= r-3; kx += 4) {
//                 float k_v0 = kf[kx     + r];
//                 float k_v1 = kf[kx + 1 + r];
//                 float k_v2 = kf[kx + 2 + r];
//                 float k_v3 = kf[kx + 3 + r];

//                 RGBAf p0 = i[y * w + (x + kx)];
//                 RGBAf p1 = i[y * w + (x + kx + 1)];
//                 RGBAf p2 = i[y * w + (x + kx + 2)];
//                 RGBAf p3 = i[y * w + (x + kx + 3)];

//                 red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
//                 green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
//                 blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
//             }

//             for (; kx <= r; kx++) {
//                 float k_v = kf[kx + r];
//                 RGBAf p = i[y * w + (x + kx)];
//                 red_val   += p.r * k_v;
//                 green_val += p.g * k_v;
//                 blue_val  += p.b * k_v;
//             }

//             RGBAf *o_p = &tf[x * h + y]; o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i[y * w + x].a;
//         }
//     }

//     for (int y = r; y < h - r; y++) {
//         for (int x = r; x < w - r; x++) {

//             float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//             int ky;

//             for (ky = -r; ky <= r-3; ky += 4) {

//                 float k_v0 = kf[ky + r];
//                 float k_v1 = kf[ky + r + 1];
//                 float k_v2 = kf[ky + r + 2];
//                 float k_v3 = kf[ky + r + 3];

//                 RGBAf p0 = tf[x * h + y + ky];
//                 RGBAf p1 = tf[x * h + y + ky + 1];
//                 RGBAf p2 = tf[x * h + y + ky + 2];
//                 RGBAf p3 = tf[x * h + y + ky + 3];

//                 red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
//                 green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
//                 blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
//             }

//             for (; ky <= r; ky++) {
//                 float k_v = k[ky + r];
//                 RGBAf p = tf[x*h + y + ky]; 
//                 red_val   += p.r * k_v;
//                 green_val += p.g * k_v;
//                 blue_val  += p.b * k_v;
//             }

//             RGBA* out = &o[y * w + x]; out->r = (unsigned char)red_val; out->g = (unsigned char)green_val; out->b = (unsigned char)blue_val; out->a = i[y * w + x].a;
//         }
//     }
// }

// void serial_t_p_t_u_omp(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o) {
    
//     #pragma omp parallel 
//     {
//         #pragma omp for
//         for (int y = 0; y < h; y++) {
//             for (int x = r; x < w - r; x++) {

//                 float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//                 int kx;

//                 #pragma omp simd
//                 for (kx = -r; kx <= r-3; kx += 4) {
//                     float k_v0 = kf[kx     + r];
//                     float k_v1 = kf[kx + 1 + r];
//                     float k_v2 = kf[kx + 2 + r];
//                     float k_v3 = kf[kx + 3 + r];

//                     RGBAf p0 = i[y * w + (x + kx)];
//                     RGBAf p1 = i[y * w + (x + kx + 1)];
//                     RGBAf p2 = i[y * w + (x + kx + 2)];
//                     RGBAf p3 = i[y * w + (x + kx + 3)];

//                     red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
//                     green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
//                     blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
//                 }

//                 for (; kx <= r; kx++) {
//                     float k_v = kf[kx + r];
//                     RGBAf p = i[y * w + (x + kx)];
//                     red_val   += p.r * k_v;
//                     green_val += p.g * k_v;
//                     blue_val  += p.b * k_v;
//                 }

//                 RGBAf *o_p = &tf[x * h + y]; o_p->r = red_val; o_p->g = green_val; o_p->b = blue_val; o_p->a = i[y * w + x].a;
//             }
//         }

//         #pragma omp for
//         for (int y = r; y < h - r; y++) {
//             for (int x = r; x < w - r; x++) {

//                 float red_val = 0.0, green_val = 0.0, blue_val = 0.0;

//                 int ky;

//                 #pragma omp simd
//                 for (ky = -r; ky <= r-3; ky += 4) {

//                     float k_v0 = kf[ky + r];
//                     float k_v1 = kf[ky + r + 1];
//                     float k_v2 = kf[ky + r + 2];
//                     float k_v3 = kf[ky + r + 3];

//                     RGBAf p0 = tf[x * h + y + ky];
//                     RGBAf p1 = tf[x * h + y + ky + 1];
//                     RGBAf p2 = tf[x * h + y + ky + 2];
//                     RGBAf p3 = tf[x * h + y + ky + 3];

//                     red_val   += p0.r * k_v0 + p1.r * k_v1 + p2.r * k_v2 + p3.r * k_v3;
//                     green_val += p0.g * k_v0 + p1.g * k_v1 + p2.g * k_v2 + p3.g * k_v3;
//                     blue_val  += p0.b * k_v0 + p1.b * k_v1 + p2.b * k_v2 + p3.b * k_v3;
//                 }

//                 for (; ky <= r; ky++) {
//                     float k_v = k[ky + r];
//                     RGBAf p = tf[x*h + y + ky]; 
//                     red_val   += p.r * k_v;
//                     green_val += p.g * k_v;
//                     blue_val  += p.b * k_v;
//                 }

//                 RGBA* out = &o[y * w + x]; out->r = (unsigned char)red_val; out->g = (unsigned char)green_val; out->b = (unsigned char)blue_val; out->a = i[y * w + x].a;
//             }
//         }
//     }
// }

// void serial_t_p_t_simd(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o) {

//     __m128 _0 = _mm_set_ps1(0);

//     for (int y = 0; y < h; y++) {
//         for (int x = r; x < w - r; x++) {

//             RGBAf_x1 acc; acc.sse = _0;
//             int kx;

//             for (kx = -r; kx <= r-3; kx += 4) {
//                 __m128 k_vv0 = k_vv[kx     + r];
//                 __m128 k_vv1 = k_vv[kx + 1 + r];
//                 __m128 k_vv2 = k_vv[kx + 2 + r];
//                 __m128 k_vv3 = k_vv[kx + 3 + r];

//                 __m128 p0 = _mm_loadu_ps((float*) &i[y * w + (x + kx)]);
//                 __m128 p1 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 1)]);
//                 __m128 p2 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 2)]);
//                 __m128 p3 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 3)]);

//                 p0 = _mm_mul_ps(p0, k_vv0);
//                 p1 = _mm_mul_ps(p1, k_vv1);
//                 p2 = _mm_mul_ps(p2, k_vv2);
//                 p3 = _mm_mul_ps(p3, k_vv3);

//                 __m128 s1 = _mm_add_ps(p0, p1);
//                 __m128 s2 = _mm_add_ps(p2, p3);

//                 acc.sse = _mm_add_ps(acc.sse, s1);
//                 acc.sse = _mm_add_ps(acc.sse, s2);
//             }

//             for (; kx <= r; kx++) {
//                 float k_v = kf[kx + r];
//                 RGBAf p = i[y * w + (x + kx)];
//                 acc.E[0] += p.r * k_v;
//                 acc.E[1] += p.g * k_v;
//                 acc.E[2] += p.b * k_v;
//             }

//             // maybe we could even make a temp of vectors, so we don't have to convert
//             RGBAf *o_p = &tf[x * h + y]; o_p->r = acc.E[0]; o_p->g = acc.E[1]; o_p->b = acc.E[2]; o_p->a = i[y * w + x].a;
//         }
//     }

//     for (int y = r; y < h - r; y++) {
//         for (int x = r; x < w - r; x++) {

//             RGBAf_x1 acc; acc.sse = _0;
//             int ky;

//             for (ky = -r; ky <= r-3; ky += 4) {
//                 __m128 k_vv0 = k_vv[ky + r    ];
//                 __m128 k_vv1 = k_vv[ky + r + 1];
//                 __m128 k_vv2 = k_vv[ky + r + 2];
//                 __m128 k_vv3 = k_vv[ky + r + 3];

//                 __m128 p0 = _mm_loadu_ps((float*) &tf[x * h + (y + ky)]);
//                 __m128 p1 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 1)]);
//                 __m128 p2 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 2)]);
//                 __m128 p3 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 3)]);

//                 p0 = _mm_mul_ps(p0, k_vv0);
//                 p1 = _mm_mul_ps(p1, k_vv1);
//                 p2 = _mm_mul_ps(p2, k_vv2);
//                 p3 = _mm_mul_ps(p3, k_vv3);

//                 __m128 s1 = _mm_add_ps(p0, p1);
//                 __m128 s2 = _mm_add_ps(p2, p3);

//                 acc.sse = _mm_add_ps(acc.sse, s1);
//                 acc.sse = _mm_add_ps(acc.sse, s2);
//             }

//             for (; ky <= r; ky++) {
//                 float k_v = kf[ky + r];
//                 RGBAf p = tf[x * h + (y + ky)];
//                 acc.E[0] += p.r * k_v;
//                 acc.E[1] += p.g * k_v;
//                 acc.E[2] += p.b * k_v;
//             }

//             RGBA* out = &o[y * w + x]; out->r = (unsigned char)acc.E[0]; out->g = (unsigned char)acc.E[1]; out->b = (unsigned char)acc.E[2]; out->a = i[y * w + x].a;
//         }
//     }
// }

// void serial_t_p_t_s_omp(int r, RGBAf* i, RGBAd* td, RGBAf* tf, int w, int h, double* k, float* kf, __m128* k_vv, int k_w, RGBA* o) {

//     __m128 _0 = _mm_set_ps1(0);

//     #pragma omp parallel 
//     {

//         #pragma omp for
//         for (int y = 0; y < h; y++) {
//             for (int x = r; x < w - r; x++) {

//                 RGBAf_x1 acc; acc.sse = _0;
//                 int kx;

//                 // #pragma omp simd
//                 for (kx = -r; kx <= r-3; kx += 4) {
//                     __m128 k_vv0 = k_vv[kx     + r];
//                     __m128 k_vv1 = k_vv[kx + 1 + r];
//                     __m128 k_vv2 = k_vv[kx + 2 + r];
//                     __m128 k_vv3 = k_vv[kx + 3 + r];

//                     __m128 p0 = _mm_loadu_ps((float*) &i[y * w + (x + kx)]);
//                     __m128 p1 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 1)]);
//                     __m128 p2 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 2)]);
//                     __m128 p3 = _mm_loadu_ps((float*) &i[y * w + (x + kx + 3)]);

//                     p0 = _mm_mul_ps(p0, k_vv0);
//                     p1 = _mm_mul_ps(p1, k_vv1);
//                     p2 = _mm_mul_ps(p2, k_vv2);
//                     p3 = _mm_mul_ps(p3, k_vv3);

//                     __m128 s1 = _mm_add_ps(p0, p1);
//                     __m128 s2 = _mm_add_ps(p2, p3);

//                     acc.sse = _mm_add_ps(acc.sse, s1);
//                     acc.sse = _mm_add_ps(acc.sse, s2);
//                 }

//                 for (; kx <= r; kx++) {
//                     float k_v = kf[kx + r];
//                     RGBAf p = i[y * w + (x + kx)];
//                     acc.E[0] += p.r * k_v;
//                     acc.E[1] += p.g * k_v;
//                     acc.E[2] += p.b * k_v;
//                 }

//                 // maybe we could even make a temp of vectors, so we don't have to convert
//                 RGBAf *o_p = &tf[x * h + y]; o_p->r = acc.E[0]; o_p->g = acc.E[1]; o_p->b = acc.E[2]; o_p->a = i[y * w + x].a;
//             }
//         }

//         #pragma omp for
//         for (int y = r; y < h - r; y++) {
//             for (int x = r; x < w - r; x++) {

//                 RGBAf_x1 acc; acc.sse = _0;
//                 int ky;

//                 // #pragma omp simd
//                 for (ky = -r; ky <= r-3; ky += 4) {
//                     __m128 k_vv0 = k_vv[ky + r    ];
//                     __m128 k_vv1 = k_vv[ky + r + 1];
//                     __m128 k_vv2 = k_vv[ky + r + 2];
//                     __m128 k_vv3 = k_vv[ky + r + 3];

//                     __m128 p0 = _mm_loadu_ps((float*) &tf[x * h + (y + ky)]);
//                     __m128 p1 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 1)]);
//                     __m128 p2 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 2)]);
//                     __m128 p3 = _mm_loadu_ps((float*) &tf[x * h + (y + ky + 3)]);

//                     p0 = _mm_mul_ps(p0, k_vv0);
//                     p1 = _mm_mul_ps(p1, k_vv1);
//                     p2 = _mm_mul_ps(p2, k_vv2);
//                     p3 = _mm_mul_ps(p3, k_vv3);

//                     __m128 s1 = _mm_add_ps(p0, p1);
//                     __m128 s2 = _mm_add_ps(p2, p3);

//                     acc.sse = _mm_add_ps(acc.sse, s1);
//                     acc.sse = _mm_add_ps(acc.sse, s2);
//                 }

//                 for (; ky <= r; ky++) {
//                     float k_v = kf[ky + r];
//                     RGBAf p = tf[x * h + (y + ky)];
//                     acc.E[0] += p.r * k_v;
//                     acc.E[1] += p.g * k_v;
//                     acc.E[2] += p.b * k_v;
//                 }

//                 RGBA* out = &o[y * w + x]; out->r = (unsigned char)acc.E[0]; out->g = (unsigned char)acc.E[1]; out->b = (unsigned char)acc.E[2]; out->a = i[y * w + x].a;
//             }
//         }
//     }
// }

// int test (Image* image) {

//     // check first row of image
//     for (int i=0; i < image->width; i += 4) {
//         RGBAf p = image->pixels[i];

//         RGBAf_x1 v; v.sse = _mm_loadu_ps((float*) &image->pixels[i]);

//         // printf("normal: %f %f %f %f\n", p.r, p.g, p.b, p.a);
//         // printf("sse   : %f %f %f %f\n", v.E[0], v.E[1], v.E[2], v.E[3]);

//         if (p.r != v.E[0] || p.g != v.E[1] || p.b != v.E[2] || p.a != v.E[3]) return 0;
//     }

//     return 1;
// }