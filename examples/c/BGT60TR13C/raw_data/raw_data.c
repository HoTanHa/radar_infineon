/* ===========================================================================
** Copyright (C) 2021-2022 Infineon Technologies AG
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are met:
**
** 1. Redistributions of source code must retain the above copyright notice,
**    this list of conditions and the following disclaimer.
** 2. Redistributions in binary form must reproduce the above copyright
**    notice, this list of conditions and the following disclaimer in the
**    documentation and/or other materials provided with the distribution.
** 3. Neither the name of the copyright holder nor the names of its
**    contributors may be used to endorse or promote products derived from
**    this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
** AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
** IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
** ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
** LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
** INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
** CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
** ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
** POSSIBILITY OF SUCH DAMAGE.
** ===========================================================================
*/

/**
 * @file    raw_data.c
 *
 * @brief   Raw data example.
 *
 * This example illustrates how to fetch time-domain data from an FMCW
 * radar sensor like BGT60TR13, BGT60UTR11AIP or BGT60ATR24 using the Radar SDK.
 */

/*
==============================================================================
   1. INCLUDE FILES
==============================================================================
*/

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "ifxAvian/DeviceControl.h"
#include "ifxBase/Base.h"

/*
==============================================================================
   2. LOCAL DEFINITIONS
==============================================================================
*/

#define NUM_FETCHED_FRAMES 2 /**< Number of frames to fetch */

/*
==============================================================================
   6. LOCAL FUNCTIONS
==============================================================================
*/
#include <math.h>
#include <complex.h>

#define PI_X ((float)(3.14159265358979323846))
const double PI = 3.14159f;

typedef double complex cplx;
typedef struct
{
    float r;
    float i;
} complex_x;

typedef struct
{
    uint32_t col_raw;
    uint32_t row_raw;
    float *data_raw;

    float *window_blackman_harris;
    float *window_chebyshev;

} my_handle_s;

static my_handle_s my_handle_data;
static my_handle_s *pHandle = &my_handle_data;
#if 0
// static void real_to_complex(float in_r[], complex_x out_c[], uint32_t len)
// {
//     for (uint32_t i = 0; i < len; i++)
//     {
//         out_c[i].r = in_r[i];
//         out_c[i].i = 0;
//     }
// }

static void init_blackmanharris(float *win, uint32_t len)
{
    if (len == 1)
    {
        win[0] = 1;
    }
    else
    {
        const float a0 = (float)0.35875;
        const float a1 = (float)0.48829;
        const float a2 = (float)0.14128;
        const float a3 = (float)0.01168;

        const float scale = PI_X / (len - 1.0f);

        for (uint32_t i = 0; i < len; ++i)
        {
            const float phi = scale * (float)i;
            win[i] = a0 + a2 * cosf(4 * phi) - a1 * cosf(2 * phi) - a3 * cosf(6 * phi);
        }
    }
}

static float acos1mx(float x)
{
    if (x < (float)0.01)
    {
        /* Use Taylor series of acos(1-x) for x=~0. The maximum relative error
         * for double is ~7e-14 at x=0.01.
         */
        const float x2 = x * x;   // x^2
        const float x3 = x * x2;  // x^3
        const float x4 = x2 * x2; // x^4
        return sqrtf(2 * x) * (1 + x / 12 + 3 * x2 / 160 + 5 * x3 / 896 + 35 * x4 / 18432);
    }
    else
        return acosf(1 - x);
}

//----------------------------------------------------------------------------

static float acosh1px(float x)
{
    return log1pf(x + sqrtf(x * (x + 2)));
}

//----------------------------------------------------------------------------

static float cosm1(float x)
{
    float s = sinf(x / 2);
    return -2 * s * s;
}

//----------------------------------------------------------------------------

static float coshm1(float x)
{
    float s = sinhf(x / 2);
    return 2 * s * s;
}

//----------------------------------------------------------------------------

static float chebyxp1(int n, float x)
{
    if (x < 0 && x >= -2)
        return cosf(n * acos1mx(-x));
    else
        return coshf(n * acosh1px(x));
}

static void init_chebyshev(float *win, uint32_t len, float at_dB)
{
    const uint32_t N = len;

    if (N == 1)
    {
        win[0] = 1.0f;
        return;
    }
    const float pssbnd_ripple = powf(10, at_dB / 20);
    const float x0m1 = coshm1(acoshf(pssbnd_ripple) / (N - 1));

    float max_val = 0;
    for (uint32_t j = 0; j < N / 2 + 1; j++)
    {
        const float n = j + (float)0.5;

        // inside the loop: sign = (-1)^i
        int sign = -1;

        float sum = 0;
        for (uint32_t i = 1; i < (N - 1) / 2 + 1; i++)
        {
            // cm1 = cos(i*pi/N)-1
            float cm1 = cosm1(i * IFX_PI / N);

            // arg = x0*cos(i*pi/N) - 1
            float arg = cm1 + x0m1 * (1 + cm1);
            sum += sign * chebyxp1(N - 1, arg) * COS((2 * IFX_PI * n * i) / N);
            sign = -sign;
        }

        // here a loss of significance occurs
        float value = pssbnd_ripple + 2 * sum;
        if (value > max_val)
            max_val = value;

        win[j] = win[N - j - 1] = value;
    }

    // ifx_vec_scale_r(win, 1 / max_val, win);
    float scale = 1 / max_val;
    for (uint32_t i = 0; i < N; ++i)
    {
        win[i] = win[i] * scale;
    }
}

void _fft(cplx buf[], cplx out[], int n, int step)
{
    if (step < n)
    {
        _fft(out, buf, n, step * 2);
        _fft(out + step, buf + step, n, step * 2);

        for (int i = 0; i < n; i += 2 * step)
        {
            cplx t = cexp(-I * PI * i / n) * out[i + step];
            buf[i / 2] = out[i] + t;
            buf[(i + n) / 2] = out[i] - t;
        }
    }
}
void fft(cplx buf[], int n)
{
    cplx *out = (cplx *)malloc(n * sizeof(cplx));
    for (int i = 0; i < n; i++)
        out[i] = buf[i];

    _fft(buf, out, n, 1);

    free(out);
}
void show(const char *s, cplx buf[], int n)
{
    printf("%s", s);
    for (int i = 0; i < n; i++)
        printf("(%10.6f %10.6f) ", creal(buf[i]), cimag(buf[i]));
    // if (!cimag(buf[i]))
    //     printf("(%10.6f  ", creal(buf[i]));
    // else
    //     printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
    printf("\r\n");
}

void print_calculate(cplx buf[], int n)
{
    int index = 32;
    double max = 0;
    double value;
    for (int i = 0; i < n; i++)
    {
        value = cabs(buf[i]) * cos(carg(buf[i]));
        printf("%10.6f  ", value);
        // printf("%8.4f ", cabs(buf[i]));

        // if (i > 32 && i < (n - 10) && (max < cabs(buf[i])))
        // {
        //     index = i;
        //     max = cabs(buf[i]);
        // }

        if (i > 32 && i < (n - 10) && (max < value))
        {
            index = i;
            max = value;
        }
    }
    printf("\r\n\n");
    printf("\tMAX  index:%d   value:%g", index, max);
    printf("\n\n");
    // for (int i = 0; i < n; i++)
    //     printf("%8.4f ", carg(buf[i]));

    // printf("\n\n");
}

float my_sum_r(float *data, uint32_t len)
{
    float sum = 0;
    for (uint32_t c = 0; c < len; c++)
    {
        sum += data[c];
    }
    return sum;
}

static void my_sub_r(float *data_in, uint32_t len, float *data_out)
{
    float sum = my_sum_r(data_in, len);
    float mean = sum / len;
    for (uint32_t c = 0; c < len; c++)
    {
        data_out[c] = data_in[c] - mean;
    }
}

static void my_mul_r(float *data_in, uint32_t len, float *window, float *data_out)
{
    for (uint32_t c = 0; c < len; c++)
    {
        data_out[c] = data_in[c] * window[c];
    }
}

static void my_fill_zero_r(float *data_in, uint32_t len_in, float *data_out, uint32_t len_out)
{
    for (uint32_t c = 0; c < len_in; c++)
    {
        data_out[c] = data_in[c];
    }
    for (uint32_t c = len_in; c < len_out; c++)
    {
        data_out[c] = 0;
    }
}

static void my_process_antenna_data(const ifx_Matrix_R_t *antenna_data)
{
    uint32_t row_raw = IFX_MAT_ROWS(antenna_data);
    uint32_t col_raw = IFX_MAT_COLS(antenna_data);

    /* get raw data */
    for (uint32_t r = 0; r < IFX_MAT_ROWS(antenna_data); r++)
    {
        for (uint32_t c = 0; c < IFX_MAT_COLS(antenna_data); c++)
        {
            ifx_Float_t value = IFX_MAT_AT(antenna_data, r, c);
            pHandle->data_raw[r * col_raw + c] = value;
        }
    }
    /* sub  -> mul -> fft fill zero -> fft*/
    float *buf_zero_pad = (float *)malloc(col_raw * 4);
    for (uint32_t r = 0; r < row_raw; r++)
    {
        float *data_sub = pHandle->data_raw + (r * col_raw);
        my_sub_r(data_sub, col_raw, data_sub);

        my_mul_r(data_sub, col_raw, pHandle->window_blackman_harris, data_sub);

        my_fill_zero_r(data_sub, col_raw, buf_zero_pad, col_raw * 4);
    }
    free(buf_zero_pad);
}
#endif

/**
 * @brief Helper function to process data
 *
 * This function is an example showing a possible way
 * of processing antenna signal. The goal in this example is to
 * sum up all chirps into one vector.
 *
 * @param antenna_data data from one antenna containing multiple chirps
 */
static void process_antenna_data(const ifx_Matrix_R_t *antenna_data)
{
    ifx_Vector_R_t chirp = {0};
    // Create the sum vector
    ifx_Vector_R_t *sum = ifx_vec_create_r(IFX_MAT_COLS(antenna_data));

    // for (uint32_t r = 0; r < IFX_MAT_ROWS(antenna_data); r++)
    // {
    //     float x = 0.0f;
    //     for (uint32_t c = 0; c < IFX_MAT_COLS(antenna_data); c++)
    //     {
    //         ifx_Float_t value = IFX_MAT_AT(antenna_data, r, c);
    //         printf("%8.4f ", value);
    //         x += value;
    //     }

    //     printf("\r\n   .........%8.4f\n", (x / IFX_MAT_COLS(antenna_data)));
    // }

    // for (uint32_t r = 0; r < IFX_MAT_ROWS(antenna_data); r++)
    // {
    //     float x = 0.0f;
    //     for (uint32_t c = 0; c < IFX_MAT_COLS(antenna_data); c++)
    //     {
    //         ifx_Float_t value = IFX_MAT_AT(antenna_data, c, r);
    //         x += value;
    //     }

    //     printf("%8.4f ", (x / IFX_MAT_COLS(antenna_data)));
    // }
    // printf("\n\n");

    // Iterate through all chirps
    for (uint32_t i = 0; i < IFX_MAT_ROWS(antenna_data); i++)
    {
        // Fetch a chirp from the antenna data matrix
        ifx_mat_get_rowview_r(antenna_data, i, &chirp);
        // add it to the sum vector
        ifx_vec_add_r(&chirp, sum, sum);
    }

    // Divide the sum vector element wise by number of chirps in the antenna data
    ifx_vec_scale_r(sum, 1.0f / IFX_MAT_ROWS(antenna_data), sum);

    uint32_t a = IFX_MAT_ROWS(antenna_data);
    uint32_t b = IFX_MAT_COLS(antenna_data);
    uint32_t asum = IFX_MAT_ROWS(sum);
    uint32_t bsum = IFX_MAT_COLS(sum);
    uint32_t achirp = IFX_MAT_ROWS(&chirp);
    uint32_t bchirp = IFX_MAT_COLS(&chirp);
    printf("process_antenna_data.....IFX_MAT_ROWS: %u  ..IFX_MAT_COLS: %u   sum:%u-%u  chirp:%u-%u\r\n",
           a, b, asum, bsum, achirp, bchirp);

    for (uint32_t i = 0; i < IFX_VEC_LEN(sum); i++)
    {
        printf("%10.6f  ", IFX_VEC_AT(sum, i));
    }
    printf("\r\n\n");

    int N = IFX_VEC_LEN(sum); // 64;

    // int len = N;              // 64;
    // float *xn1 = (float *)malloc(N * sizeof(float));
    // float *Xr = (float *)malloc(N * sizeof(float));
    // float *Xi = (float *)malloc(N * sizeof(float));
    // float *ampli = (float *)malloc(N * sizeof(float));
    // uint32_t *ampli_pp = (uint32_t *)malloc(N * sizeof(uint32_t));

    // for (uint32_t i = 0; i < IFX_VEC_LEN(sum); i++)
    // {
    //     xn1[i] = IFX_VEC_AT(sum, i);
    // }
    // for (int k = 0; k < N; k++)
    // {
    //     Xr[k] = 0;
    //     Xi[k] = 0;
    //     for (int n = 0; n < len; n++)
    //     {
    //         Xr[k] = (Xr[k] + xn1[n] * cos(2 * 3.14159265 * k * n / N));
    //         Xi[k] = (Xi[k] - xn1[n] * sin(2 * 3.14159265 * k * n / N));
    //     }
    //     ampli[k] = sqrt(Xr[k] * Xr[k] + Xi[k] * Xi[k]);
    //     ampli_pp[k] = (uint32_t)(1000 * ampli[k]);
    //     printf("%8.4f ", ampli[k]);
    //     // printf("aaaa=%g  (%g) + j(%g)\n", xn1[k], Xr[k], Xi[k]);
    // }
    // printf("\n");
    // uint32_t dinh = ampli_pp[9];
    // int index = 9;
    // for (int k = 9; k < (N - 1); k++)
    // {
    //     if (dinh < ampli_pp[k])
    //     {
    //         dinh = ampli_pp[k];
    //         index = k;
    //     }
    // }
    // printf("\tMax value: %u    ...index: %d\n\n", dinh, index);

    // free(xn1);
    // free(Xr);
    // free(Xi);
    // free(ampli);
    // free(ampli_pp);

    int sizefft = N; // 4 * N;
    cplx *Xfft = (cplx *)malloc(sizefft * sizeof(cplx));

    for (uint32_t i = 0; i < (uint32_t)(sizefft); i++)
    {
        Xfft[i] = 0.0f;
    }
    for (uint32_t i = 0; i < IFX_VEC_LEN(sum); i++)
    {
        Xfft[i] = IFX_VEC_AT(sum, i);
        // Xfft[4 * i + 1] = IFX_VEC_AT(sum, i);
        // Xfft[4 * i + 2] = IFX_VEC_AT(sum, i);
        // Xfft[4 * i + 3] = IFX_VEC_AT(sum, i);
    }
    // // show("Data: \n", (cplx *)Xfft, sizefft);
    // fft((cplx *)Xfft, sizefft);
    // // fft((cplx *)Xfft, sizefft);
    // // show("FFT : \n", (cplx *)Xfft, sizefft);

    // print_calculate(Xfft, sizefft);
    free(Xfft);

    printf("\n\n");

    ifx_vec_destroy_r(sum);
}

//----------------------------------------------------------------------------

/**
 * @brief Helper function to separate signals
 *
 * Separates different antenna signals
 * and pass them for further processing.
 *
 * @param frame The frame may contain multiple antenna signals,
 *              depending on the device configuration.
 *              Each antenna signal can contain multiple chirps.
 */
static void process_frame(ifx_Cube_R_t *frame)
{
    ifx_Matrix_R_t antenna_data;

    static int count = 0;
    count++;
    uint32_t a = IFX_CUBE_ROWS(frame);
    printf("%5d........process_frame...............: %u \r\n", count, a);

    // for (uint32_t i = 0; i < a; i++)
    // {
    //     ifx_cube_get_row_r(frame, i, &antenna_data);
    //     process_antenna_data(&antenna_data);
    // }

    ifx_cube_get_row_r(frame, 0, &antenna_data);
    process_antenna_data(&antenna_data);
}

/**
 * @brief Helper function to print device configuration
 *
 * Print the device configuration to stdout.
 *
 * @param config        device configuration
 */
static void print_device_config(const ifx_Avian_Config_t *config)
{
    printf("Device configuration:\n");
    printf("sample_rate_Hz:          %" PRIu32 "\n", config->sample_rate_Hz);
    printf("rx_mask:                 %" PRIu32 "\n", config->rx_mask);
    printf("tx_mask:                 %" PRIu32 "\n", config->tx_mask);
    printf("tx_power_level:          %" PRIu32 "\n", config->tx_power_level);
    printf("if_gain_dB:              %" PRIu32 "\n", config->if_gain_dB);
    printf("start_frequency_Hz:      %" PRIu64 "\n", config->start_frequency_Hz);
    printf("end_frequency_Hz:        %" PRIu64 "\n", config->end_frequency_Hz);
    printf("num_samples_per_chirp:   %" PRIu32 "\n", config->num_samples_per_chirp);
    printf("num_chirps_per_frame:    %" PRIu32 "\n", config->num_chirps_per_frame);
    printf("chirp_repetition_time_s: %g\n", config->chirp_repetition_time_s);
    printf("frame_repetition_time_s: %g\n", config->frame_repetition_time_s);
    printf("hp_cutoff_Hz:            %" PRIu32 "\n", config->hp_cutoff_Hz);
    printf("aaf_cutoff_Hz:           %" PRIu32 "\n", config->aaf_cutoff_Hz);
    printf("mimo_mode:               %s\n", config->mimo_mode == IFX_MIMO_TDM ? "time-domain multiplexed" : "off");
    printf("\n");
}

/*
==============================================================================
   7. MAIN METHOD
==============================================================================
 */

int main(int argc, char **argv)
{
    ifx_Error_t error = IFX_OK;
    ifx_Avian_Config_t device_config = {0};
    ifx_Avian_Device_t *device_handle = NULL;
    ifx_Cube_R_t *frame = NULL;

    printf("Radar SDK Version: %s\n", ifx_sdk_get_version_string_full());

    /* Open the device: Connect to the first radar sensor found. */
    device_handle = ifx_avian_create();
    if ((error = ifx_error_get()) != IFX_OK)
    {
        fprintf(stderr, "Failed to open device: %s\n", ifx_error_to_string(error));
        goto out;
    }

    const char *uuid = ifx_avian_get_board_uuid(device_handle);
    printf("UUID of board: %s\n", uuid);

    /* Get default device configuration for connected radar sensor.
     * Here, we use the default configuration that depends on the connected
     * radar sensor.
     * Typically, the radar sensor is either configured by setting all members
     * of the device configuration (ifx_Avian_Config_t) or by setting the
     * members of a structure of type ifx_Avian_Metrics_t and converting the
     * structure to a device configuration using the function
     * ifx_avian_metrics_to_config. The latter approach using the
     * ifx_Avian_Metrics_t structure provides an easy way to configure the
     * radar sensor in terms of high-level parameters, while the device
     * configuration structure ifx_Avian_Config_t allows to configure the
     * sensor in a more detailed way.
     */
    ifx_avian_get_config_defaults(device_handle, &device_config);
    if ((error = ifx_error_get()) != IFX_OK)
    {
        fprintf(stderr, "Failed to get default device config:  %s\n", ifx_error_to_string(error));
        goto out;
    }

    /**************************************************************/
    // device_config.num_samples_per_chirp = 128;
    // device_config.num_chirps_per_frame = 128;
    // device_config.sample_rate_Hz = 2000000;
    // device_config.end_frequency_Hz = 61000000000;
    // device_config.start_frequency_Hz = 60000000000;
    /*********************************************************************/

    /* Apply the device settings based on the device configuration structure. */
    ifx_avian_set_config(device_handle, &device_config);
    if ((error = ifx_error_get()) != IFX_OK)
    {
        fprintf(stderr, "Failed to set device config:  %s\n", ifx_error_to_string(error));
        goto out;
    }

    /* Print the current device configuration to stdout:
     * Read back the current device configuration. Due to discrete register
     * values in the radar sensors, the device configuration passed to
     * ifx_avian_set_config and the actual device configuration set might
     * differ slightly.
     * We then print the read-back device configuration to stdout.
     */
    ifx_avian_get_config(device_handle, &device_config);
    print_device_config(&device_config);

    uint64_t deltaF = (device_config.end_frequency_Hz - device_config.start_frequency_Hz) / 1000000;
    float deltaR = 3 * 100 / (2.0f * deltaF);
    float Tc = device_config.num_samples_per_chirp * 1000000.0f / device_config.sample_rate_Hz;
    float Knum = deltaF * 1.0f / (device_config.num_samples_per_chirp * 1.0f / device_config.sample_rate_Hz);
    // float Knum = deltaF * 1.0f / (device_config.chirp_repetition_time_s);
    float Rmax = (device_config.sample_rate_Hz * 3 * 100) / (2 * Knum);

    printf("Bandwidth = BW = %ld*10^6Hz \r\n", deltaF);
    printf("deltaR = c/(2*BW) = %.4f \r\n", deltaR);
    printf("Chirp duration = Tc = %.4f * 10^-6 (second)\r\n", Tc);
    printf("Slope K = BW/Tc = %.4f * 10^6 \r\n", Knum);
    printf("Rmax = Fs*c/2Knum = %.4f ==> [-%4fm %4fm]\r\n", Rmax, Rmax / 2, Rmax / 2);
    float frame_rate = 1 / device_config.frame_repetition_time_s;
    printf("frame_rate = 1/frame_repetition_time = %.4f \r\n", frame_rate);
    uint64_t fCenter = (device_config.end_frequency_Hz / 1000000 + device_config.start_frequency_Hz / 1000000) / 2;
    float wavelength = 3 * 100.0f / fCenter;
    printf("f_center = %ld *10^6 (Hz)   ==> wavelength = c/f_center = %.4f (m) \r\n", fCenter, wavelength);
    float r_res = deltaR;
    float r_max = r_res * device_config.num_samples_per_chirp / 2;
    printf("r_max = r_res * dc.num_samples_per_chirp / 2 =  %.4f(m)\r\n", r_max);
    float v_max = wavelength / (4 * device_config.chirp_repetition_time_s);
    printf("v_max = wavelength / (4 * dc.chirp_repetition_time_s) = %.4f(m/s)\r\n", v_max);
    float v_res = 2 * v_max / device_config.num_chirps_per_frame;
    printf("v_res = 2 * v_max / dc.num_chirps_per_frame = %.4f(m/s)\r\n", v_res);

    printf("\r\n\r\n");
    sleep(5);

    pHandle->col_raw = device_config.num_samples_per_chirp;
    pHandle->row_raw = device_config.num_chirps_per_frame;
    pHandle->data_raw = (float *)malloc(pHandle->col_raw * pHandle->row_raw * sizeof(float));

    pHandle->window_blackman_harris = (float *)malloc(device_config.num_samples_per_chirp * sizeof(float));
    pHandle->window_chebyshev = (float *)malloc(device_config.num_chirps_per_frame * sizeof(float));

    // init_blackmanharris(pHandle->window_blackman_harris, device_config.num_samples_per_chirp);
    // init_chebyshev(pHandle->window_chebyshev, device_config.num_chirps_per_frame, 100);

    /* Fetch NUM_FETCHED_FRAMES number of frames. */
    for (int frame_number = 0; frame_number < NUM_FETCHED_FRAMES; frame_number++)
    {
        /* Get the time-domain data for the next frame. The function will block
         * until the full frame is available and copy the data into the frame
         * handle.
         * This function also creates a frame structure for time domain data
         * acquisition, if not created already. It is the responsibility of
         * the caller to free the returned frame in this scope.
         */
        frame = ifx_avian_get_next_frame(device_handle, frame);
        if ((error = ifx_error_get()) != IFX_OK)
        {
            fprintf(stderr, "Failed to get next frame: %s\n", ifx_error_to_string(error));
            goto out;
        }

        /* Process the frame. */
        process_frame(frame);
        usleep(1000);
    }

    free(pHandle->data_raw);
    free(pHandle->window_blackman_harris);
    free(pHandle->window_chebyshev);

out:
    /* Close the device after processing all frames. It is valid to pass NULL
     * to destroy functions.
     */
    ifx_avian_destroy(device_handle);
    ifx_cube_destroy_r(frame);

    return error == IFX_OK ? EXIT_SUCCESS : EXIT_FAILURE;
}
