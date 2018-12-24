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

#define OUTPUT_FILENAME "sound_level_meter.m"
#define INPUT_FILENAME "sound.dat"
#define INPUT_FILENAME_LIMIT 100

// Constants used in rest of the program
#define OCT_ULIMIT_CENTER_FREQ   0.45f
#define OCT_HALFBAND_FACTOR      1.1225f /** 2^1/6 */
#define OCT_HALF_PI              M_PI * 0.5
#define OCT_EDGE_FREQS           30
#define OCT_PIVOT_FREQ           1.0f/20.0f
#define OCT_ERROR_CODE           -1
#define OCT_CUBEROOT_TWO         OCT_HALFBAND_FACTOR * OCT_HALFBAND_FACTOR /** 2^1/3 */
#define OCT_CUBEROOT_TWO_SQ      OCT_CUBEROOT_TWO * OCT_CUBEROOT_TWO /** 2^2/3 */

// Typical sampling rates are supported but 8KHz would be too
// low for a meaningful analysis
#define DEFAULT_SAMPLING_RATE    44100
#define MIN_REQ_SAMPLING_RATE    11025
// 30 seconds of Audio at the default sampling rate is the maximum size supported.
#define MAX_INPUT_DATA_SIZE   DEFAULT_SAMPLING_RATE * 30

// print usage/help message
void usage()
{
    printf("octave3filt -- 1/3 octave filter design\n");
    printf("options (default values in []):\n");
    printf("  u/h   : print usage/help\n");
    printf("  n     : filter order, n > 0 [3]\n");
    printf("  r     : sampling rate, must be higher than 11025Hz, [44.1KHz]\n");
    printf("  d     : filename containing 32 bit FP sound data sampled at r [sound.dat] \n");
}

// data_f: pointer to array of 32bit FP, containing the data
// filename: File from which to read floating point data
// returns number of samples that were read
size_t get_sound_data(float* data_f, char* filename)
{
    uint64_t samples_read = OCT_ERROR_CODE;
    if (filename) {
        FILE* fid = fopen(filename, "r");
        if (fid != NULL) {
            samples_read = 0;
            while(samples_read <= MAX_INPUT_DATA_SIZE) {
                if (fscanf(fid, "%f", &data_f[samples_read]) == EOF) {
                    fprintf(stderr, "Info: Finished parsing input file %s with samples %zu\n",
                            filename, samples_read);
                    break;
                } else {
                    samples_read++;
                }
            }
            fclose(fid);
        } else {
            fprintf(stderr, "Failed to read file %s\n", filename);
        }
    } else {
        fprintf(stderr, "Invalid filename/path supplied\n");
    }

    return samples_read;
}
// n: Filter order
// f0: Center frequency
// b: Numerator coefficients of digital IIR filter, length 2*N+1
// a: Denominator coefficient of digital IIR filter, length 2*N+1
// returns -1 on Error
int one_third_filter(unsigned int n, float f0, float* b, float* a) {

    // validate input
    if (f0 <= 0 || f0 >= OCT_ULIMIT_CENTER_FREQ) {
        fprintf(stderr,"error: %12.4e, centre frequency out of range\n", f0);
        usage();
        return OCT_ERROR_CODE;
    }

    if (b == NULL || a == NULL) {
        fprintf(stderr, "Invalid co-efficient's array supplied to one-third octave design function\n");
        return OCT_ERROR_CODE;
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
    float f1 = f0 / OCT_HALFBAND_FACTOR;
    float f2 = f0 * OCT_HALFBAND_FACTOR;
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
    return 0;
}

int main(int argc, char*argv[]) {

    // Default IIR filter order
    unsigned int n = 3;
    // Default sampling rate
    float fs = DEFAULT_SAMPLING_RATE;
    // Default input data filename
    char input_file[INPUT_FILENAME_LIMIT];
    strncpy(input_file, INPUT_FILENAME, strlen(INPUT_FILENAME)+1);

    // Boiler plate argument parsing
    int dopt;
    while ((dopt = getopt(argc,argv,"uh:n:r:c:d:")) != EOF) {
        switch (dopt) {
        case 'u':
        case 'h':
            usage();
            return 0;
        case 'r':
            fs = atof(optarg);
            if (fs < MIN_REQ_SAMPLING_RATE) {
                fprintf(stderr, "error: Too low sampling rate for a range of sound analysis");
                exit(1);
            }
            break;
        case 'n':
            n  = atoi(optarg);
            if (n == 0 || n > 10) {
                fprintf(stderr, "error: Invalid filter order supplied \n");
                exit(1);
            }
            break;
        case 'd':
            if (strlen(optarg) > INPUT_FILENAME_LIMIT) {
                fprintf(stderr, "Error: Too many characters in input file name, limit 100 characters\n");
                usage();
                exit(1);
            }
            strncpy(input_file, optarg, strlen(optarg)+1);
            fprintf(stdout, "Input file supplied %s\n", input_file);
            break;
        default:
            fprintf(stderr, "No valid option supplied\n");
            usage();
            exit(1);
        }
    }
    // Construct frequency set, while keeping highest upper edge freq band
    // sufficiently lower than sampling rate.
    float band_edges[OCT_EDGE_FREQS] = {25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630,
                        800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000,
                        12500, 16000, 20000};
    float* freq_set = (float *) malloc(OCT_EDGE_FREQS * sizeof(float));
    unsigned f_index = 0;
    int e_index;
    for (e_index = (-OCT_EDGE_FREQS/2) - 1; e_index < (OCT_EDGE_FREQS/2) - 1; e_index++) {
        float center_freq = 1000.0f * pow(OCT_HALFBAND_FACTOR * OCT_HALFBAND_FACTOR, e_index);
        if (center_freq >= fs/3) {
            fprintf(stderr,"Warning!, Octave analysis will be restricted till %12.4e=\n",
                    band_edges[f_index]);
            break;
        }
        freq_set[f_index] = center_freq / fs;
        fprintf(stdout, "freq-edge: %d= normalized_edge: %12.4e=\n", f_index, freq_set[f_index]);
        f_index++;
    }
    assert(f_index < OCT_EDGE_FREQS);
    fprintf(stdout, "Info: Maximum band edge that will be processed %12.4e=\n", band_edges[f_index]);

    // Find pivot frequency, below which octave filters will be applied
    // after decimation.
    unsigned p_index;
    for (p_index = 0; p_index < f_index; p_index++) {
        if (freq_set[p_index] > OCT_PIVOT_FREQ) {
            fprintf(stdout, "Info, Pivot index=%d and associated frequency %12.4e=\n",
                    p_index, freq_set[p_index]);
            break;
        }
        p_index++;
    }

    // Read sound data from input file
    float* sound_data = (float *) malloc(MAX_INPUT_DATA_SIZE * sizeof(float));
    size_t data_read = get_sound_data(sound_data, input_file);
    if (data_read == -1) {
        fprintf(stderr, "Error reading from file %s \n", input_file);
        free(freq_set);
        free(sound_data);
        exit(1);
    }

    // We have the data and the centre frequenices, now we can start 1/3 octave
    // filterbank analysis
    float b[2*n + 1];
    float a[2*n + 1];
    // RMS power of each band, initialized to NAN
    float* rms_power = (float *) malloc(OCT_EDGE_FREQS * sizeof(float));
    for (unsigned rr = 0; rr < OCT_EDGE_FREQS; rr++) {
        rms_power[rr] = NAN;
    }

    // Above pivot frequency, we can simply apply 1/3 octave filter in each band
    fprintf(stdout,"Info: Applying 1/3 Octave filter band from freq_set[%d]=%12.8e"
            "till freq_set[%d]=%12.8e\n", p_index, freq_set[p_index], f_index, freq_set[f_index]);

    for (unsigned b_index = p_index; b_index <= f_index; b_index++)
    {
        if (one_third_filter(n, freq_set[b_index], b, a) == OCT_ERROR_CODE) {
            fprintf(stderr, "Failed to design octave 1/3 filter for freq_set[%d]=%12.8e\n",
                    b_index, freq_set[b_index]);
            free(freq_set);
            free(rms_power);
            free(sound_data);
            exit(1);
        }
        iirfilt_rrrf octave_one_third = iirfilt_rrrf_create(b, 2*n + 1, a, 2*n + 1);
        float accum_power = 0.0f;
        float filt_output = 0.0f;
        for (size_t data_index = 0; data_index < data_read; data_index++) {
            iirfilt_rrrf_execute(octave_one_third, sound_data[data_index], &filt_output);
            accum_power += filt_output * filt_output;
        }
        rms_power[b_index] = 10*log10f(accum_power);
        fprintf(stdout, "Info: Power computed rms_power[%d]=%12.8f\n", b_index, rms_power[b_index]);
    }

    // For lower frequencies we apply same filter-band with three frequencies
    // albeit with decimation.
    float bu[2*n + 1];
    float au[2*n + 1];
    float bc[2*n + 1];
    float ac[2*n + 1];
    float bl[2*n + 1];
    float al[2*n + 1];
    float half_fc = freq_set[p_index-1]/2;
    if (one_third_filter(n, half_fc, bu, au) == OCT_ERROR_CODE) {
        fprintf(stderr, "Failed to design octave 1/3 filter for freq=%12.8e\n", half_fc);
        free(freq_set);
        free(rms_power);
        free(sound_data);
        exit(1);
    }
    if (one_third_filter(n, half_fc/OCT_CUBEROOT_TWO, bc, ac) == OCT_ERROR_CODE) {
        fprintf(stderr, "Failed to design octave 1/3 filter for freq=%12.8e\n", half_fc/OCT_CUBEROOT_TWO);
        free(freq_set);
        free(rms_power);
        free(sound_data);
        exit(1);
    }
    if (one_third_filter(n, half_fc/OCT_CUBEROOT_TWO_SQ, bl, al) == OCT_ERROR_CODE) {
        fprintf(stderr, "Failed to design octave 1/3 filter for freq=%12.8e\n", half_fc/OCT_CUBEROOT_TWO_SQ);
        free(freq_set);
        free(rms_power);
        free(sound_data);
        exit(1);
    }

    free(freq_set);
    free(rms_power);
    free(sound_data);
    return 0;
}
