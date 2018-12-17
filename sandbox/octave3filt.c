// octave3filt.c
//
// Tests 1/3 Octave filter design.
// Based on IIR butter-worth filter example
// Original work by TBA

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include "liquid.h"

#define OUTPUT_FILENAME "octave3filt.m"


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


int main(int argc, char*argv[]) {
    // options
    unsigned int n=3;       // filter order
    float f0 = 8000.0f/44100.0f;       // center frequency (band-pass, band-stop)

    int dopt;
    while ((dopt = getopt(argc,argv,"uht:b:n:r:s:f:c:o:")) != EOF) {
        switch (dopt) {
        case 'u':
        case 'h':
            usage();
            return 0;
        case 'n': n = atoi(optarg);         break;
        case 'c': f0 = atof(optarg);        break;
        default:
            exit(1);
        }
    }

    // validate input
    if (f0 <= 0 || f0 >= OCT_ULIMIT_CENTER_FREQ) {
        fprintf(stderr,"error: %s, centre frequency out of range\n", argv[0]);
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

    unsigned int i;

    printf("Butterworth filter design:\n");
    nza = 0;
    k0 = 1.0f;
    butter_azpkf(n,za,pa,&ka);

    printf("poles (analog):\n");
    for (i=0; i<npa; i++)
        printf("  pa[%3u] = %12.8f + j*%12.8f\n", i, crealf(pa[i]), cimagf(pa[i]));
    printf("zeros (analog):\n");
    for (i=0; i<nza; i++)
        printf("  za[%3u] = %12.8f + j*%12.8f\n", i, crealf(za[i]), cimagf(za[i]));
    printf("gain (analog):\n");
    printf("  ka : %12.8f + j*%12.8f\n", crealf(ka), cimagf(ka));

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
    printf("  alpha = %12.4e, Qd = %12.4e, Qr = %12.4e, f2 = %12.4e, f1 = %12.4e;\n", alpha, Qd, Qr, f2, f1);
    assert(alpha > 1.0f);
    float fc = f0 * alpha;
    float m = iirdes_freqprewarp(LIQUID_IIRDES_BANDPASS,fc,f0);
    printf("m : %12.8f\n", m);
    bilinear_zpkf(za,    nza,
                  pa,    npa,
                  k0,    m,
                  zd, pd, &kd);

    // open output file
    FILE*fid = fopen(OUTPUT_FILENAME,"w");
    fprintf(fid,"%% %s : auto-generated file\n", OUTPUT_FILENAME);
    fprintf(fid,"clear all;\n");
    fprintf(fid,"close all;\n");

    printf("zeros (digital, low-pass prototype):\n");
    for (i=0; i<n; i++)
        printf("  zd[%3u] = %12.4e + j*%12.4e;\n", i, crealf(zd[i]), cimagf(zd[i]));
    printf("poles (digital, low-pass prototype):\n");
    for (i=0; i<n; i++)
        printf("  pd[%3u] = %12.4e + j*%12.4e;\n", i, crealf(pd[i]), cimagf(pd[i]));
    printf("gain (digital):\n");
    printf("  kd : %12.8f + j*%12.8f\n", crealf(kd), cimagf(kd));

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

    // update paramteres : n -> 2*n
    n = 2*n;

    fprintf(fid,"f0=%12.4e\n",f0);
    fprintf(fid,"n=%u;\n", n);
    fprintf(fid,"nfft=1024;\n");

    // print digital z/p/k
    fprintf(fid,"zd = zeros(1,n);\n");
    fprintf(fid,"pd = zeros(1,n);\n");
    for (i=0; i<n; i++) {
        fprintf(fid,"  zd(%3u) = %12.4e + j*%12.4e;\n", i+1, crealf(zd[i]), cimagf(zd[i]));
        fprintf(fid,"  pd(%3u) = %12.4e + j*%12.4e;\n", i+1, crealf(pd[i]), cimagf(pd[i]));
    }


    float b[n+1];       // numerator
    float a[n+1];       // denominator

    // convert complex digital poles/zeros/gain into transfer function
    iirdes_dzpk2tff(zd,pd,n,kd,b,a);

    // print coefficients
    for (i=0; i<=n; i++) printf("a[%3u] = %12.8f;\n", i, a[i]);
    for (i=0; i<=n; i++) printf("b[%3u] = %12.8f;\n", i, b[i]);

    fprintf(fid,"a = zeros(1,n+1);\n");
    fprintf(fid,"b = zeros(1,n+1);\n");
    for (i=0; i<=n; i++) {
        fprintf(fid,"a(%3u) = %12.4e;\n", i+1, a[i]);
        fprintf(fid,"b(%3u) = %12.4e;\n", i+1, b[i]);
    }
    fprintf(fid,"\n");
    fprintf(fid,"H = fft(b,nfft)./fft(a,nfft);\n");
    fprintf(fid,"H = fftshift(H);\n");
    fprintf(fid,"%% group delay\n");
    fprintf(fid,"c = conv(b,fliplr(conj(a)));\n");
    fprintf(fid,"cr = c.*[0:(length(c)-1)];\n");
    fprintf(fid,"t0 = fftshift(fft(cr,nfft));\n");
    fprintf(fid,"t1 = fftshift(fft(c, nfft));\n");
    fprintf(fid,"polebins = find(abs(t1)<1e-6);\n");
    fprintf(fid,"t0(polebins)=0;\n");
    fprintf(fid,"t1(polebins)=1;\n");
    fprintf(fid,"gd = real(t0./t1) - length(a) + 1;\n");

    // plot zeros, poles
    fprintf(fid,"\n");
    fprintf(fid,"figure;\n");
    fprintf(fid,"k=0:0.01:1;\n");
    fprintf(fid,"ti = cos(2*pi*k);\n");
    fprintf(fid,"tq = sin(2*pi*k);\n");
    fprintf(fid,"plot(ti,tq,'-','LineWidth',1,'Color',[1 1 1]*0.7,...\n");
    fprintf(fid,"     real(zd),imag(zd),'o','LineWidth',2,'Color',[0.5 0   0],'MarkerSize',2,...\n");
    fprintf(fid,"     real(pd),imag(pd),'x','LineWidth',2,'Color',[0   0.5 0],'MarkerSize',2);\n");
    fprintf(fid,"xlabel('real');\n");
    fprintf(fid,"ylabel('imag');\n");
    fprintf(fid,"title('z-plane');\n");
    fprintf(fid,"grid on;\n");
    fprintf(fid,"axis([-1 1 -1 1]*1.2);\n");
    fprintf(fid,"axis square;\n");

    // plot group delay
    fprintf(fid,"f = [0:(nfft-1)]/nfft - 0.5;\n");
    fprintf(fid,"figure;\n");
    fprintf(fid,"  plot(f,gd,'-','Color',[0 0.5 0],'LineWidth',2);\n");
    fprintf(fid,"  axis([0.0 0.5 0 ceil(1.1*max(gd))]);\n");
    fprintf(fid,"  grid on;\n");
    fprintf(fid,"  xlabel('Normalized Frequency');\n");
    fprintf(fid,"  ylabel('Group delay [samples]');\n");

    // plot magnitude response
    fprintf(fid,"figure;\n");
    fprintf(fid,"subplot(2,1,1),\n");
    fprintf(fid,"  plot(f,20*log10(abs(H)),'-','Color',[0.5 0 0],'LineWidth',2);\n");
    fprintf(fid,"  axis([0.0 0.5 -4 1]);\n");
    fprintf(fid,"  grid on;\n");
    fprintf(fid,"  xlabel('Normalized Frequency');\n");
    fprintf(fid,"  ylabel('Filter PSD [dB]');\n");
    fprintf(fid,"subplot(2,1,2),\n");
    fprintf(fid,"  plot(f,20*log10(abs(H)),'-','Color',[0.5 0 0],'LineWidth',2);\n");
    fprintf(fid,"  axis([0.0 0.5 -100 10]);\n");
    fprintf(fid,"  grid on;\n");
    fprintf(fid,"  xlabel('Normalized Frequency');\n");
    fprintf(fid,"  ylabel('Filter PSD [dB]');\n");

    fclose(fid);
    printf("results written to %s.\n", OUTPUT_FILENAME);

    printf("done.\n");
    return 0;
}

