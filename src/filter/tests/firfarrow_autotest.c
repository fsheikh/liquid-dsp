/*
 * Copyright (c) 2007 - 2015 Joseph Gaeddert
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "autotest/autotest.h"
#include "liquid.h"

//
// AUTOTEST : fir farrow filter, 3rd order, n=8 co-efficients
//
void autotest_firfarrow_q3_n8()
{
    // create input params
    unsigned int n = 8;
    unsigned int q = 3;
    float fc = 3150.0f / 9600.0f;
    float As = 60.0f;

    float tol = 1e-3f;

    firfarrow_rrrf farrow = firfarrow_rrrf_create(n, q, fc, As);
    firfarrow_rrrf_set_delay(farrow, 0.0f);
    firfarrow_rrrf_print(farrow);
    // create testing vectors
    float exp_impulse_res[8] = {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    float cal_impulse_res[8];
    firfarrow_rrrf_get_coefficients(farrow, &cal_impulse_res[0]);

    unsigned int i = 0;
    for (i=0; i<n; i++) {
        CONTEND_DELTA( cal_impulse_res[i], exp_impulse_res[i], tol );
    }

    // destroy filter
    firfarrow_rrrf_destroy(farrow);
}


//
// AUTOTEST : iir group delay, n=3
//
#if 0
void autotest_iir_groupdelay_n3()
{
    // create coefficients array
    float b[3] = {0.20657210,  0.41314420, 0.20657210};
    float a[3] = {1.00000000, -0.36952737, 0.19581573};

    float tol = 1e-3f;
    unsigned int i;

    // create testing vectors
    float fc[4] = { 0.000,
                    0.125,
                    0.250,
                    0.375};
    
    float g0[4] = { 0.973248939389634,
                    1.366481121240365,
                    1.227756735863196,
                    0.651058521306726};

    // run tests
    float g;
    for (i=0; i<4; i++) {
        g = iir_group_delay(b, 3, a, 3, fc[i]);
        CONTEND_DELTA( g, g0[i], tol );
    }

    // create filter
    iirfilt_rrrf filter = iirfilt_rrrf_create(b,3,a,3);

    // run tests again
    for (i=0; i<4; i++) {
        g = iirfilt_rrrf_groupdelay(filter, fc[i]);
        CONTEND_DELTA( g, g0[i], tol );
    }

    // destroy filter
    iirfilt_rrrf_destroy(filter);
}

#endif
