// sound_level_meter.c
//
// Software implementation of a sound level meter
// Based on original work in Matlab/Octave by 
//

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
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
#define OCT_ERROR_CODE           -1
#define OCT_CUBEROOT_TWO         (OCT_HALFBAND_FACTOR * OCT_HALFBAND_FACTOR) /** 2^1/3 */
#define OCT_CUBEROOT_TWO_SQ      (OCT_CUBEROOT_TWO * OCT_CUBEROOT_TWO) /** 2^2/3 */
#define OCT_FIRFILT_DELAY        7
#define OCT_DECIMATION_FACTOR    2
// Typical sampling rates are supported but 8KHz would be too
// low for a meaningful analysis
#define DEFAULT_SAMPLING_RATE    44100
#define MIN_REQ_SAMPLING_RATE    11025
// 30 seconds of Audio at the default sampling rate is the maximum size supported.
#define MAX_INPUT_DATA_SIZE     DEFAULT_SAMPLING_RATE * 30
#define OCT_BLK_SIZE            64
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

// data_f: pointer to array of 32bit FP samples of sound signal
// filename: File from which to read floating point sound data
// returns number of samples that were read
size_t get_sound_data(float* data_f, char* filename)
{
    uint64_t samples_read = OCT_ERROR_CODE;
    if (filename) {
        FILE* fid = fopen(filename, "r");
        if (fid != NULL) {
            samples_read = 0;
            while(samples_read <= MAX_INPUT_DATA_SIZE) {
                if (fscanf(fid, "%e", &data_f[samples_read]) == EOF) {
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
// fs: Sampling rate for which filter is designed.
// b: Numerator coefficients of digital IIR filter, length 2*N+1
// a: Denominator coefficient of digital IIR filter, length 2*N+1
// returns -1 on Error
int one_third_filter(unsigned int n, float f0, float fs, float* b, float* a) {

    // validate input
    if (f0 <= 0 || f0 >= OCT_ULIMIT_CENTER_FREQ * fs) {
        fprintf(stderr,"error: %12.4e, centre frequency out of range\n", f0);
        return OCT_ERROR_CODE;
    }

    if (b == NULL || a == NULL) {
        fprintf(stderr, "Invalid co-efficient's array supplied to one-third octave design function\n");
        return OCT_ERROR_CODE;
    }

    // complex digital poles/zeros/gain
    // NOTE: allocated double the filter order to cover band-pass, band-stop cases
    float f1 = f0 / OCT_HALFBAND_FACTOR;
    float f2 = f0 * OCT_HALFBAND_FACTOR;
    float Qr = f0 / (f2 - f1);
    float Qd = Qr * (OCT_HALF_PI / n) / sinf(OCT_HALF_PI / n);
    float alpha = ((1 + sqrtf(1 + (4 * powf(Qd, 2.0f)))) / 2.0f) / Qd;
    assert(alpha > 1.0f);
    printf("alpha=%12.4e fc=%12.4e f0=%12.4e f2=%12.4e\r\n", alpha, f0*alpha, f0, f2);
    liquid_iirdes(LIQUID_IIRDES_BUTTER,
                  LIQUID_IIRDES_BANDPASS,
                  LIQUID_IIRDES_TF,
                  n,
                  alpha*(f0/fs), f0/fs, 60.0f, 1.0f,
                  b, a);
    return 0;
}

unsigned construct_freq_set(float* freq_set, unsigned* pivot, float sampling_rate)
{
    if (freq_set == NULL) {
        fprintf(stderr, "Invalid frequency set vector received\n");
        return OCT_ERROR_CODE;
    }
    // Construct frequency set, while keeping highest upper edge freq band
    // sufficiently lower than sampling rate.
    float band_edges[OCT_EDGE_FREQS] = {25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630,
                        800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000,
                        12500, 16000, 20000};
    // Conversion from band-edges to band-pass center frequency
    unsigned f_index = 0;
    for (int e_index = (-OCT_EDGE_FREQS/2) - 1; e_index < (OCT_EDGE_FREQS/2) - 1; e_index++) {
        freq_set[f_index] = 1000.0f * powf(1.25999f, e_index);
        f_index++;
    }
    if (sampling_rate < 3 * band_edges[OCT_EDGE_FREQS - 1]) {
        fprintf(stderr, "Info! sampling rate not enough for processing whole frequeny range\n");
        f_index = 0;
        while (freq_set[f_index] < sampling_rate/3.0f) {
            freq_set[f_index] = freq_set[f_index];
            f_index++ ;
        }
        if (f_index < OCT_EDGE_FREQS - 1) {
            fprintf(stderr,"Warning!, Octave analysis will be restricted till band_edge[%d]=%12.4e\n",
                    f_index, band_edges[f_index]);
        }
    }
    assert(f_index < OCT_EDGE_FREQS);
    // Find pivot frequency, below which octave filters will be applied
    // after decimation.
    unsigned p_index = 0;
    while ((p_index < f_index) && (freq_set[p_index] < sampling_rate/20.0f)) {
        p_index++;
    }
    printf("final index %d, pivot index %d\n", f_index, p_index);
    *pivot = p_index;
    return f_index;
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
                fprintf(stderr, "error: Too low sampling rate for a range of sound analysis\n");
                exit(1);
            }
            break;
        case 'n':
            n  = atoi(optarg);
            if (n == 0 || n > 10) {
                fprintf(stderr, "error: Invalid filter order supplied\n");
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

    // Read sound data from input file
    float* sound_data = (float *) malloc(MAX_INPUT_DATA_SIZE * sizeof(float));
    size_t data_read = get_sound_data(sound_data, input_file);
    if (data_read == -1) {
        fprintf(stderr, "Error reading from file %s \n", input_file);
        free(sound_data);
        exit(1);
    }

    // Output file to verify with GNU OCTAVE
    FILE* fid = fopen(OUTPUT_FILENAME, "w");
    // RMS power of each band, initialized to NAN
    float rms_power[OCT_EDGE_FREQS];
    for (unsigned rr = 0; rr < OCT_EDGE_FREQS; rr++) {
        rms_power[rr] = NAN;
    }
    fprintf(fid, "%% %s : auto-generated file\n", OUTPUT_FILENAME);
    fprintf(fid, "close all;\n");
    fprintf(fid, "clear all;\n");
    fprintf(fid, "rms_power=NaN*ones(1,%d);\n", OCT_EDGE_FREQS);
    fprintf(fid, "freq_set=zeros(1,%d);\n", OCT_EDGE_FREQS);
    fprintf(fid, "num_samples=%lu;\n", data_read);
    fprintf(fid, "XX=[%12.4e %12.4e %12.4e %12.4e %12.4e];\n", sound_data[9], sound_data[99], sound_data[999], sound_data[1999], sound_data[2999]);
    fprintf(fid, "Fs=%12.4e;\n", fs);

    float freq_set[OCT_EDGE_FREQS]; // Set of center frequencies derived from 1/3 octave band edges
    unsigned p_index; // pivot frequency index
    unsigned f_index = construct_freq_set(freq_set, &p_index, fs);
    for (unsigned fii=0; fii < OCT_EDGE_FREQS; fii++) {
        if (fii == p_index) {
            fprintf(fid, "pivot_freq=%12.8e;\n", freq_set[fii] * fs);
        }
        if (fii == f_index) {
            fprintf(fid, "final_freq=%12.8e;\n", freq_set[fii] * fs);
        }
        fprintf(fid, "freq_set(%d)=%12.8e;\n", fii+1, freq_set[fii] * fs);
    }

    // We have the data and the centre frequencies, now we can start 1/3 octave
    // filterbank analysis
    unsigned L = 2*n + 1;
    fprintf(fid, "L=%d\n", L);
    for (unsigned b_index = p_index; b_index < f_index; b_index++)
    {
        float b[L];
        float a[L];
        fprintf(fid, "b_%d=zeros(1, %d);\n", b_index+1, L);
        fprintf(fid, "a_%d=zeros(1, %d);\n", b_index+1, L);

        if (one_third_filter(n, freq_set[b_index], fs, b, a) == OCT_ERROR_CODE) {
            fprintf(stderr, "Failed to design octave 1/3 filter for freq_set[%d]=%12.8e\n",
                    b_index, freq_set[b_index]);
            free(sound_data);
            exit(1);
        }
        iirfilt_rrrf octave_one_third = iirfilt_rrrf_create(b, L, a, L);
        for (unsigned findex = 0; findex < L; findex++) {
            fprintf(fid, "b_%d(%d)=%12.4e;\n", b_index+1, findex+1, b[findex]);
            fprintf(fid, "a_%d(%d)=%12.4e;\n", b_index+1, findex+1, a[findex]);
        }
#if 0
        float accum_power = 0.0f;
        float filt_output[OCT_BLK_SIZE] = {0.0f};
        for (size_t data_index = 0; data_index < data_read; data_index+=OCT_BLK_SIZE) {
            iirfilt_rrrf_execute_block(octave_one_third, &sound_data[data_index], OCT_BLK_SIZE, filt_output);
            accum_power += liquid_sumsqf(filt_output, OCT_BLK_SIZE);
        }
#endif
        float accum_power = 0.0f;
        float filt_output = 0.0f;
        for (size_t data_index = 0; data_index < data_read; data_index++) {
            iirfilt_rrrf_execute(octave_one_third, sound_data[data_index], &filt_output);
            accum_power += liquid_sumsqf(&filt_output, 1);
        }
        rms_power[b_index] = 10*log10f(accum_power/(float) data_read);
        fprintf(fid, "rms_power(%d)=%12.8f\n", b_index+1, rms_power[b_index]);
        iirfilt_rrrf_destroy(octave_one_third);
    }
    // Call GNU octave for verification.
    fprintf(fid, "x=load('%s');\n",input_file);
    fprintf(fid, "[P,F]=filtbank(x, %12.4e, %12.4e, '%s');\n",fs, data_read/fs, "extended"); 
    // For lower frequencies we apply same filter-band with three frequencies
    // albeit with decimation.
    float bu[2*n + 1];
    float au[2*n + 1];
    float bc[2*n + 1];
    float ac[2*n + 1];
    float bl[2*n + 1];
    float al[2*n + 1];
    float half_fc = freq_set[p_index-1] * 2.0f;
    if (one_third_filter(n, half_fc, fs, bu, au) == OCT_ERROR_CODE) {
        fprintf(stderr, "Failed to design upper octave 1/3 filter for freq=%12.8e\n", half_fc);
        free(sound_data);
        fclose(fid);
        exit(1);
    }
    if (one_third_filter(n, half_fc/OCT_CUBEROOT_TWO, fs, bc, ac) == OCT_ERROR_CODE) {
        fprintf(stderr, "Failed to design center octave 1/3 filter for freq=%12.8e\n", half_fc/OCT_CUBEROOT_TWO);
        free(sound_data);
        exit(1);
    }
    fprintf(stdout, "upper freq=%12.4e, center_freq=%12.4e, lower_freq=%12.4e\n", half_fc, half_fc/OCT_CUBEROOT_TWO, half_fc/OCT_CUBEROOT_TWO_SQ);
    if (one_third_filter(n, half_fc/OCT_CUBEROOT_TWO_SQ, fs, bl, al) == OCT_ERROR_CODE) {
        fprintf(stderr, "Failed to design lower octave 1/3 filter for freq=%12.8e\n", half_fc/OCT_CUBEROOT_TWO_SQ);
        free(sound_data);
        fclose(fid);
        exit(1);
    }
    // FIR filter used for decimation
    firdecim_rrrf decim = firdecim_rrrf_create_kaiser(OCT_DECIMATION_FACTOR, OCT_FIRFILT_DELAY, 30.0f);

    firdecim_rrrf_set_scale(decim, 1.0/2.0f);

    size_t decimated_size = data_read;
    float filt_output[OCT_BLK_SIZE] = {0.0f};
    // Lower bands are treated with same set of a three-set filterbank
    for (int b_index = p_index - 1; b_index >= 0; b_index = b_index-3) {
        // Create decimated output
        decimated_size = decimated_size / 2;
        float* decimated_out = (float *) malloc(decimated_size * sizeof(float));
        firdecim_rrrf_execute_block(decim, sound_data, decimated_size, decimated_out);
        sound_data = (float *) realloc(sound_data, decimated_size * sizeof(float));
        if (sound_data == NULL) {
            fprintf(stderr, "Cannot resize input after decimation\n");
            free(decimated_out);
            firdecim_rrrf_destroy(decim);
            fclose(fid);
            exit(1);
        }
        memcpy(sound_data, decimated_out, decimated_size * sizeof(float));
        printf("Resized input data to %lu samples\n", decimated_size);
        for (size_t di=0; di < decimated_size; di++) {
            if (isnan(sound_data[di])) {
                fprintf(stderr, "Info: NaN found in sound data after decimation at index %lu\n", di);
                free(decimated_out);
                firdecim_rrrf_destroy(decim);
                free(decimated_out);
                fclose(fid);
                exit(1);
            }
        }
        free(decimated_out);
        // Apply octave filters analysis on decimated output
        iirfilt_rrrf octave_one_third_upper = iirfilt_rrrf_create(bu, 2*n + 1, au, 2*n + 1);
        iirfilt_rrrf octave_one_third_center = iirfilt_rrrf_create(bc, 2*n + 1, ac, 2*n + 1);
        iirfilt_rrrf octave_one_third_lower = iirfilt_rrrf_create(bl, 2*n + 1, al, 2*n + 1);
        float accum_power[3]  = {0.0f, 0.0f, 0.0f};
        bool apply_center = b_index >= 1 ? true : false;
        bool apply_lower = b_index >= 2 ? true : false;
        for (size_t ii=0; ii < decimated_size; ii+=OCT_BLK_SIZE) {
            // Start with upper ocatve filter
            iirfilt_rrrf_execute_block(octave_one_third_upper, &sound_data[ii], OCT_BLK_SIZE, filt_output);
            if(isnan(liquid_sumsqf(filt_output, OCT_BLK_SIZE))) {
                fprintf(stdout, "NAN after applying upper 1/3 octave filter for sound_data from index=%lu\n", ii);
            }
            accum_power[2] += liquid_sumsqf(filt_output, OCT_BLK_SIZE);
            iirfilt_rrrf_reset(octave_one_third_upper);
            for (size_t fi=0; fi < OCT_BLK_SIZE; fi++) {
                filt_output[fi] = 0.0f;
            }
            // Continue to center octave filter if possible
            if (apply_center) {
                iirfilt_rrrf_execute_block(octave_one_third_center, &sound_data[ii], OCT_BLK_SIZE, filt_output);
                accum_power[1] += liquid_sumsqf(filt_output, OCT_BLK_SIZE);
                if(isnan(accum_power[1])) {
                    fprintf(stdout, "NAN after applying center 1/3 octave filter for sound_data from index=%lu\n", ii);
                }
                iirfilt_rrrf_reset(octave_one_third_center);
                for (size_t fi=0; fi < OCT_BLK_SIZE; fi++) {
                    filt_output[fi] = 0.0f;
                }
            }
            // Finally attempt lowest octave filter if possible
            if (apply_lower) {
                iirfilt_rrrf_execute_block(octave_one_third_lower, &sound_data[ii], OCT_BLK_SIZE, filt_output);
                accum_power[0] += liquid_sumsqf(filt_output, OCT_BLK_SIZE);
                if(isnan(accum_power[0])) {
                    fprintf(stdout, "NAN after applying lower 1/3 octave filter for sound_data from index=%lu\n", ii);
                }
                iirfilt_rrrf_reset(octave_one_third_lower);
                for (size_t fi=0; fi < OCT_BLK_SIZE; fi++) {
                    filt_output[fi] = 0.0f;
                }
            }
            //printf("Dying in which outer loop %d b_index and inner index %zu\n", b_index, ii);
        }
        rms_power[b_index] = 10*log10f(accum_power[2]/(float) decimated_size);
        accum_power[2] = 0.0f;
        fprintf(stdout, "Info: Power computed rms_power[%d]=%12.8f from %lu samples\n", b_index, rms_power[b_index], decimated_size);
        if (apply_center) {
            rms_power[b_index - 1] = 10*log10f(accum_power[1]/(float) decimated_size);
            accum_power[1] = 0.0f;
            fprintf(stdout, "Info: Power computed rms_power[%d]=%12.8f from %lu samples\n", b_index - 1, rms_power[b_index - 1], decimated_size);
        }
        if (apply_lower) {
            rms_power[b_index - 2] =10*log10f(accum_power[0]/(float) decimated_size);
            accum_power[0] = 0.0f;
            fprintf(stdout, "Info: Power computed rms_power[%d]=%12.8f from %lu samples\n", b_index - 2, rms_power[b_index-2], decimated_size);
        }
        fprintf(stdout, "Info: Center frequency processed: freq_set[%d]=%12.4e\n", b_index, freq_set[b_index]);
        iirfilt_rrrf_destroy(octave_one_third_upper);
        iirfilt_rrrf_destroy(octave_one_third_center);
        iirfilt_rrrf_destroy(octave_one_third_lower);
    }
    // Add the followin the output GNU octave file
    //>> axis([0 31 0 40])
    //>> set(gca, 'XTickLabel', F)
    //>> bar(rms_power)
    // Cleanup owned resources by the program
    firdecim_rrrf_destroy(decim);
    fprintf(stdout, "Debug: Liquid resources cleared\n");
    free(sound_data);
    fclose(fid);
    return 0;
}
