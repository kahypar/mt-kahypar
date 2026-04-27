#include <stdlib.h>
#include <immintrin.h>
#include <math.h>

#undef __cplusplus

#include "mt-kahypar/partition/coarsening/multilevel/ml/feature_definitions.h"


static inline void kernel_32_48(const float *A, const float *B, const float *bias, float *C)
{
    __m256 c00 = _mm256_load_ps(bias + 0);
    __m256 c01 = _mm256_load_ps(bias + 8);
    __m256 c02 = _mm256_load_ps(bias + 16);
    __m256 c03 = _mm256_load_ps(bias + 24);
    __m256 c04 = _mm256_load_ps(bias + 32);
    __m256 c05 = _mm256_load_ps(bias + 40);
    __m256 c10 = _mm256_load_ps(bias + 0);
    __m256 c11 = _mm256_load_ps(bias + 8);
    __m256 c12 = _mm256_load_ps(bias + 16);
    __m256 c13 = _mm256_load_ps(bias + 24);
    __m256 c14 = _mm256_load_ps(bias + 32);
    __m256 c15 = _mm256_load_ps(bias + 40);

    for (int i = 0; i < 32; i++)
    {
        __m256 b0 = _mm256_load_ps(B + i * 48 + 0);
        __m256 b1 = _mm256_load_ps(B + i * 48 + 8);
        __m256 b2 = _mm256_load_ps(B + i * 48 + 16);
        __m256 b3 = _mm256_load_ps(B + i * 48 + 24);
        __m256 b4 = _mm256_load_ps(B + i * 48 + 32);
        __m256 b5 = _mm256_load_ps(B + i * 48 + 40);

        __m256 a0 = _mm256_broadcast_ss(A + (32 * 0) + i);
        c00 = _mm256_fmadd_ps(b0, a0, c00);
        c01 = _mm256_fmadd_ps(b1, a0, c01);
        c02 = _mm256_fmadd_ps(b2, a0, c02);
        c03 = _mm256_fmadd_ps(b3, a0, c03);
        c04 = _mm256_fmadd_ps(b4, a0, c04);
        c05 = _mm256_fmadd_ps(b5, a0, c05);

        __m256 a1 = _mm256_broadcast_ss(A + (32 * 1) + i);
        c10 = _mm256_fmadd_ps(b0, a1, c10);
        c11 = _mm256_fmadd_ps(b1, a1, c11);
        c12 = _mm256_fmadd_ps(b2, a1, c12);
        c13 = _mm256_fmadd_ps(b3, a1, c13);
        c14 = _mm256_fmadd_ps(b4, a1, c14);
        c15 = _mm256_fmadd_ps(b5, a1, c15);
    }

    c00 = _mm256_max_ps(c00, _mm256_setzero_ps());
    c01 = _mm256_max_ps(c01, _mm256_setzero_ps());
    c02 = _mm256_max_ps(c02, _mm256_setzero_ps());
    c03 = _mm256_max_ps(c03, _mm256_setzero_ps());
    c04 = _mm256_max_ps(c04, _mm256_setzero_ps());
    c05 = _mm256_max_ps(c05, _mm256_setzero_ps());
    c10 = _mm256_max_ps(c10, _mm256_setzero_ps());
    c11 = _mm256_max_ps(c11, _mm256_setzero_ps());
    c12 = _mm256_max_ps(c12, _mm256_setzero_ps());
    c13 = _mm256_max_ps(c13, _mm256_setzero_ps());
    c14 = _mm256_max_ps(c14, _mm256_setzero_ps());
    c15 = _mm256_max_ps(c15, _mm256_setzero_ps());

    _mm256_store_ps(C + (48 * 0) + 0, c00);
    _mm256_store_ps(C + (48 * 0) + 8, c01);
    _mm256_store_ps(C + (48 * 0) + 16, c02);
    _mm256_store_ps(C + (48 * 0) + 24, c03);
    _mm256_store_ps(C + (48 * 0) + 32, c04);
    _mm256_store_ps(C + (48 * 0) + 40, c05);
    _mm256_store_ps(C + (48 * 1) + 0, c10);
    _mm256_store_ps(C + (48 * 1) + 8, c11);
    _mm256_store_ps(C + (48 * 1) + 16, c12);
    _mm256_store_ps(C + (48 * 1) + 24, c13);
    _mm256_store_ps(C + (48 * 1) + 32, c14);
    _mm256_store_ps(C + (48 * 1) + 40, c15);
}

static inline void kernel_48_32(const float *A, const float *B, const float *bias, float *C)
{
    __m256 c00 = _mm256_load_ps(bias + 0);
    __m256 c01 = _mm256_load_ps(bias + 8);
    __m256 c02 = _mm256_load_ps(bias + 16);
    __m256 c03 = _mm256_load_ps(bias + 24);
    __m256 c10 = _mm256_load_ps(bias + 0);
    __m256 c11 = _mm256_load_ps(bias + 8);
    __m256 c12 = _mm256_load_ps(bias + 16);
    __m256 c13 = _mm256_load_ps(bias + 24);

    for (int i = 0; i < 48; i++)
    {
        __m256 b0 = _mm256_load_ps(B + i * 32 + 0);
        __m256 b1 = _mm256_load_ps(B + i * 32 + 8);
        __m256 b2 = _mm256_load_ps(B + i * 32 + 16);
        __m256 b3 = _mm256_load_ps(B + i * 32 + 24);

        __m256 a0 = _mm256_broadcast_ss(A + (48 * 0) + i);
        c00 = _mm256_fmadd_ps(b0, a0, c00);
        c01 = _mm256_fmadd_ps(b1, a0, c01);
        c02 = _mm256_fmadd_ps(b2, a0, c02);
        c03 = _mm256_fmadd_ps(b3, a0, c03);

        __m256 a1 = _mm256_broadcast_ss(A + (48 * 1) + i);
        c10 = _mm256_fmadd_ps(b0, a1, c10);
        c11 = _mm256_fmadd_ps(b1, a1, c11);
        c12 = _mm256_fmadd_ps(b2, a1, c12);
        c13 = _mm256_fmadd_ps(b3, a1, c13);
    }

    c00 = _mm256_max_ps(c00, _mm256_setzero_ps());
    c01 = _mm256_max_ps(c01, _mm256_setzero_ps());
    c02 = _mm256_max_ps(c02, _mm256_setzero_ps());
    c03 = _mm256_max_ps(c03, _mm256_setzero_ps());
    c10 = _mm256_max_ps(c10, _mm256_setzero_ps());
    c11 = _mm256_max_ps(c11, _mm256_setzero_ps());
    c12 = _mm256_max_ps(c12, _mm256_setzero_ps());
    c13 = _mm256_max_ps(c13, _mm256_setzero_ps());

    _mm256_store_ps(C + (32 * 0) + 0, c00);
    _mm256_store_ps(C + (32 * 0) + 8, c01);
    _mm256_store_ps(C + (32 * 0) + 16, c02);
    _mm256_store_ps(C + (32 * 0) + 24, c03);
    _mm256_store_ps(C + (32 * 1) + 0, c10);
    _mm256_store_ps(C + (32 * 1) + 8, c11);
    _mm256_store_ps(C + (32 * 1) + 16, c12);
    _mm256_store_ps(C + (32 * 1) + 24, c13);
}

static inline void kernel_32_16(const float *A, const float *B, const float *bias, float *C)
{
    __m256 c00 = _mm256_load_ps(bias + 0);
    __m256 c01 = _mm256_load_ps(bias + 8);
    __m256 c10 = _mm256_load_ps(bias + 0);
    __m256 c11 = _mm256_load_ps(bias + 8);
    __m256 c20 = _mm256_load_ps(bias + 0);
    __m256 c21 = _mm256_load_ps(bias + 8);
    __m256 c30 = _mm256_load_ps(bias + 0);
    __m256 c31 = _mm256_load_ps(bias + 8);

    for (int i = 0; i < 32; i++)
    {
        __m256 b0 = _mm256_load_ps(B + i * 16 + 0);
        __m256 b1 = _mm256_load_ps(B + i * 16 + 8);

        __m256 a0 = _mm256_broadcast_ss(A + (32 * 0) + i);
        c00 = _mm256_fmadd_ps(b0, a0, c00);
        c01 = _mm256_fmadd_ps(b1, a0, c01);

        __m256 a1 = _mm256_broadcast_ss(A + (32 * 1) + i);
        c10 = _mm256_fmadd_ps(b0, a1, c10);
        c11 = _mm256_fmadd_ps(b1, a1, c11);

        __m256 a2 = _mm256_broadcast_ss(A + (32 * 2) + i);
        c20 = _mm256_fmadd_ps(b0, a2, c20);
        c21 = _mm256_fmadd_ps(b1, a2, c21);

        __m256 a3 = _mm256_broadcast_ss(A + (32 * 3) + i);
        c30 = _mm256_fmadd_ps(b0, a3, c30);
        c31 = _mm256_fmadd_ps(b1, a3, c31);
    }

    c00 = _mm256_max_ps(c00, _mm256_setzero_ps());
    c01 = _mm256_max_ps(c01, _mm256_setzero_ps());
    c10 = _mm256_max_ps(c10, _mm256_setzero_ps());
    c11 = _mm256_max_ps(c11, _mm256_setzero_ps());
    c20 = _mm256_max_ps(c20, _mm256_setzero_ps());
    c21 = _mm256_max_ps(c21, _mm256_setzero_ps());
    c30 = _mm256_max_ps(c30, _mm256_setzero_ps());
    c31 = _mm256_max_ps(c31, _mm256_setzero_ps());

    _mm256_store_ps(C + (16 * 0) + 0, c00);
    _mm256_store_ps(C + (16 * 0) + 8, c01);
    _mm256_store_ps(C + (16 * 1) + 0, c10);
    _mm256_store_ps(C + (16 * 1) + 8, c11);
    _mm256_store_ps(C + (16 * 2) + 0, c20);
    _mm256_store_ps(C + (16 * 2) + 8, c21);
    _mm256_store_ps(C + (16 * 3) + 0, c30);
    _mm256_store_ps(C + (16 * 3) + 8, c31);
}

static inline void kernel_16_1(const float *A, const float *B, const float *bias, float *C, int n)
{
    for (int i = 0; i < n; i++)
    {
        float res = *bias;
        for (int j = 0; j < 16; j++)
            res += A[i * 16 + j] * B[j];

        C[i] = 1.0f / (1.0f + expf(-res));
    }
}

void column_to_row_matrix(const float *in, float *out, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            out[i * m + j] = in[j * n + i];
        }
    }
}

float* precompute_params() {
    int offsets_W[] = {0, 1584, 3152, 3680};

    float *params_row = (float*)aligned_alloc(32, sizeof(float) * 3697);
    for (int i = 0; i < 3697; i++)
        params_row[i] = params[i];

    column_to_row_matrix(params + offsets_W[0], params_row + offsets_W[0], 32, 48);
    column_to_row_matrix(params + offsets_W[1], params_row + offsets_W[1], 48, 32);
    column_to_row_matrix(params + offsets_W[2], params_row + offsets_W[2], 32, 16);

    return params_row;
}

void free_params(float *params_row) {
    free(params_row);
}

void predict(float *in, float *params_row, int n, float *out)
{
    int offsets_W[] = {0, 1584, 3152, 3680};
    int offsets_B[] = {1536, 3120, 3664, 3696};

    for (int i = 0; i < n; i += 2)
    {
        kernel_32_48(in + i * 32, params_row + offsets_W[0], params_row + offsets_B[0], out + i * 48);
    }
    for (int i = 0; i < n; i += 2)
    {
        kernel_48_32(out + i * 48, params_row + offsets_W[1], params_row + offsets_B[1], in + i * 32);
    }
    for (int i = 0; i < n; i += 4)
    {
        kernel_32_16(in + i * 32, params_row + offsets_W[2], params_row + offsets_B[2], out + i * 16);
    }
    kernel_16_1(out, params_row + offsets_W[3], params_row + offsets_B[3], out, n);
}
