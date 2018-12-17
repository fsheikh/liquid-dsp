// sound_level_meter.c
//
// Software implementation of a sound level meter
// Based on original work in Matlab/Octave by 
//

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include "liquid.h"

#define OUTPUT_FILENAME "SoundLevelMeter.m"


// Constants used in rest of the program
#define OCT_ULIMIT_CENTER_FREQ   0.45
#define OCT_HALF_BAND_FACTOR     1.1225 /** 2^1/6 */
#define OCT_HALF_PI              M_PI * 0.5
// print usage/help message
void usage()
{
    printf("octave3filt -- 1/3 octave filter design\n");
    printf("options (default values in []):\n");
    printf("  u/h   : print usage/help\n");
    printf("  n     : filter order, n > 0 [3]\n");
    printf("  c     : center frequency, 0 < c < 0.5 [0.25]\n");
}

// n: Filter order
// f0: Center frequency
// b: Numerator coefficients of digital IIR filter, length 2*N+1
// a: Denominator coefficient of digital IIR filter, length 2*N+1
void one_third_filter(unsigned int n, float f0, float* b, float* a) {

    // validate input
    if (f0 <= 0 || f0 >= OCT_ULIMIT_CENTER_FREQ) {
        fprintf(stderr,"error: %12.4e, centre frequency out of range\n", f0);
        usage();
        exit(1);
    }

    // number of analaog poles/zeros
    unsigned int npa = n;
    unsigned int nza;

    // analog poles/zeros/gain
    float complex pa[n];
    float complex za[n];
    float complex ka;
    float complex k0;

    nza = 0;
    k0 = 1.0f;
    butter_azpkf(n,za,pa,&ka);

    // complex digital poles/zeros/gain
    // NOTE: allocated double the filter order to cover band-pass, band-stop cases
    float complex zd[2*n];
    float complex pd[2*n];
    float complex kd;
    float f1 = f0 / OCT_HALF_BAND_FACTOR;
    float f2 = f0 * OCT_HALF_BAND_FACTOR;
    float Qr = f0 / (f2 - f1);
    float Qd = Qr * (OCT_HALF_PI / n) / sin(OCT_HALF_PI / n);
    float alpha = (1 + sqrt(1 + 4 * powf(Qd, 2))) / 2.0f / Qd;
    assert(alpha > 1.0f);
    float fc = f0 * alpha;
    float m = iirdes_freqprewarp(LIQUID_IIRDES_BANDPASS,fc,f0);
    printf("m : %12.8f\n", m);
    bilinear_zpkf(za,    nza,
                  pa,    npa,
                  k0,    m,
                  zd, pd, &kd);

    // transform zeros, poles in band-pass since octave 1/3
    // should have bandpass characteristics
    // allocate memory for transformed zeros, poles
    float complex zd1[2*n];
    float complex pd1[2*n];

    // run zeros, poles trasform
    iirdes_dzpk_lp2bp(zd, pd,   // low-pass prototype zeros, poles
                      n,        // filter order
                      f0,       // center frequency
                      zd1,pd1); // transformed zeros, poles (length: 2*n)

    // copy transformed zeros, poles
    memmove(zd, zd1, 2*n*sizeof(float complex));
    memmove(pd, pd1, 2*n*sizeof(float complex));

    // convert complex digital poles/zeros/gain into transfer function
    iirdes_dzpk2tff(zd,pd,2*n,kd,b,a);

    // print coefficients
    unsigned int i;
    for (i=0; i<=2*n; i++) printf("a[%3u] = %12.8f;\n", i, a[i]);
    for (i=0; i<=2*n; i++) printf("b[%3u] = %12.8f;\n", i, b[i]);
}

int main(int argc, char*argv[]) {
    unsigned int n = 3;
    float b[2*n + 1];
    float a[2*n + 1];
    one_third_filter(n, 8000.0/44100.0, b, a);
    return 0;
}
