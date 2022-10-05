#include "fft.h"
#include <stdio.h>

/**
 * This file contains the implementation of `fft_compute` and ancillary
 * functions.
 */

// Include for math functions and definition of PI
#include <math.h>
// Included to get access to `malloc` and `free`
#include <stdlib.h>

// Forward declaration of helper methods
void get_even(const complex* in, complex* out, const int n)
{
    for (int i = 0; i < n / 2; i++) {
        // Transfer all the even indexed numbers to the output array
        out[i] = in[2 * i];
    }
}

void get_odd(const complex* in, complex* out, const int n)
{
    for (int i = 0; i < n / 2; i++) {
        // Transfer all the odd indexed numbers to the output array
        out[i] = in[2 * i + 1];
    }
}

// reverser en bitstreng
unsigned int rev(unsigned int n, int s){
    unsigned int r = 0;
    // s er lengden på bit-strngen
    while(s>>=1){
        r <<= 1;
        r|= n & 1;
        n >>= 1;
    }
    return r;
}   

// Iterativ Cooley Tukey
void bit_rev_copy(const complex* in, complex* out, int len)
{
    for (int i = 0; i < len; i++)
        out[rev(i, len)] = in[i];
}

void fft_compute(const complex* in, complex* out, const int n)
{
    bit_rev_copy(in, out, n);
    int m = 2; // 2^i
    int m2 = 1; // 2^(i-1)
    complex w[n>>1]; // twiddle faktorene
    w[0] = 1;
    for (int i = 1; i <= __builtin_ctz(n); i++, m<<=1, m2<<=1) {
        w[1] = cexp(-M_PI / m2 * I);
        for (int j = 2; j< m2; j++)
            w[j] = w[1]*w[j-1];
        for (int k = 0; k < n; k+=m) {
            // Denne loopen er lik den i naive versjonen
            for (int j = 0; j < m2; j++) {
                complex u = out[k + j];
                complex t = w[j] * out[k + j + m2];
                out[k + j] = u + t;
                out[k + j + m2] = u - t;
            }
        }
    }
}

// Naiv implementation:
void fft_compute_naiv(const complex* in, complex* out, const int n)
{
    if (n == 1) {
        out[0] = in[0];
    } else {
        const int half = n / 2;
        // First we declare and allocate arrays
        // Allocate enough room for half the input values
        complex* even = malloc(sizeof(complex) * half);
        complex* odd = malloc(sizeof(complex) * half);
        complex* even_out = malloc(sizeof(complex) * half);
        complex* odd_out = malloc(sizeof(complex) * half);
        // Extract even and odd indexed numbers using methods above
        get_even(in, even, n);
        get_odd(in, odd, n);
        // Recursively calculate the result for bottom and top half
        fft_compute(even, even_out, n / 2);
        fft_compute(odd, odd_out, n / 2);
        // Combine the output of the two previous recursions
        for (int i = 0; i < half; ++i) {
            const complex e = even_out[i];
            const complex o = odd_out[i];
            const complex w = cexp(0 - (2. * M_PI * i) / n * I);
            out[i] = e + w * o;
            out[i + half] = e - w * o;
        }
        // Since we allocated room for variables we need to release
        // the memory!
        free(even);
        free(odd);
        free(even_out);
        free(odd_out);
    }
}
