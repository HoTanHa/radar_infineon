

#include <string.h>
#include <unistd.h>
#include <signal.h>

#include "common.h"
#include "ifxAlgo/Algo.h"
#include "ifxAvian/Avian.h"
#include "ifxRadar/RangeAngleImage.h"

#define IFX_ADC_SAMPLERATE_HZ (1000000U)
#define IFX_RX_MASK (1)
#define IFX_TX_MASK (1)
#define IFX_BGT_TX_POWER (31U)
#define IFX_IF_GAIN_DB (33)
#define IFX_RANGE_RESOLUTION_M (0.03f)
#define IFX_MAX_RANGE_M (0.96f)
#define IFX_SPEED_RESOLUTION_M_S (0.067f)
#define IFX_SPEED_RESOLUTION_M_S1 (0.08f)
#define IFX_MAX_SPEED_M_S (2.096f)
#define IFX_FRAME_RATE_HZ (12.946f)
#define IFX_MAXIMUM_SPEED_MPS (2.5f)
#define IFX_ALPHA_MTI_FILTER (0.5f)
#define IFX_SPECT_THRESHOLD (1e-6f)
typedef struct
{
    ifx_RAI_t *rai_handle;
    ifx_Cube_R_t *rai_cube_r;
    ifx_Matrix_R_t *rdm_r;

    ifx_OSCFAR_t *oscfar_handle;
    ifx_DBSCAN_t *dbscan_handle;
} app_handle_t;

static ifx_Avian_Metrics_t default_metrics =
    {
        .range_resolution_m = IFX_RANGE_RESOLUTION_M,
        .max_range_m = IFX_MAX_RANGE_M,
        .speed_resolution_m_s = IFX_SPEED_RESOLUTION_M_S,
        .max_speed_m_s = IFX_MAX_SPEED_M_S,
        .center_frequency_Hz = 60.5e9f, // use default value for Avian
};

void rdm_peak_search(const ifx_Matrix_R_t *matrix, uint32_t *rmax, uint32_t *cmax)
{
    *rmax = 0, *cmax = 0;
    ifx_Float_t max_value = IFX_MAT_AT(matrix, 0, 0);

    for (uint32_t r = 0; r < IFX_MAT_ROWS(matrix); r++)
    {
        for (uint32_t c = 0; c < IFX_MAT_COLS(matrix); c++)
        {
            ifx_Float_t value = IFX_MAT_AT(matrix, r, c);
            if (value > max_value)
            {
                *rmax = r;
                *cmax = c;
                max_value = value;
            }
        }
    }
}

ifx_Error_t app_radar_init(app_handle_t *app_context)
{
    memset(app_context, 0, sizeof(app_handle_t));
    return ifx_error_get();
}

ifx_Error_t app_radar_config(app_handle_t *app_context, ifx_Avian_Device_t *device, ifx_json_t *json, ifx_Avian_Config_t *dev_config)
{
    ifx_Error_t ret = 0;

    const uint32_t range_fft_size = dev_config->num_samples_per_chirp * 4;  // Zero padding of 4 gives good range resolution in Range spectrum
    const uint32_t doppler_fft_size = dev_config->num_chirps_per_frame * 4; // Zero padding of 4 gives good range resolution in Doppler spectrum

    printf("range_fft_size:%u   doppler_fft_size:%u\r\n", range_fft_size, doppler_fft_size);

    printf("Device configuration:\n");
    printf("sample_rate_Hz:          %u\n", dev_config->sample_rate_Hz);
    printf("rx_mask:                 %u\n", dev_config->rx_mask);
    printf("tx_mask:                 %u\n", dev_config->tx_mask);
    printf("tx_power_level:          %u\n", dev_config->tx_power_level);
    printf("if_gain_dB:              %u\n", dev_config->if_gain_dB);
    printf("start_frequency_Hz:      %lu\n", dev_config->start_frequency_Hz);
    printf("end_frequency_Hz:        %lu\n", dev_config->end_frequency_Hz);
    printf("num_samples_per_chirp:   %u\n", dev_config->num_samples_per_chirp);
    printf("num_chirps_per_frame:    %u\n", dev_config->num_chirps_per_frame);
    printf("chirp_repetition_time_s: %g\n", dev_config->chirp_repetition_time_s);
    printf("frame_repetition_time_s: %g\n", dev_config->frame_repetition_time_s);
    printf("hp_cutoff_Hz:            %u\n", dev_config->hp_cutoff_Hz);
    printf("aaf_cutoff_Hz:           %u\n", dev_config->aaf_cutoff_Hz);
    printf("\n");

    ifx_PPFFT_Config_t range_fft_config = (ifx_PPFFT_Config_t){
        .fft_type = IFX_FFT_TYPE_R2C,
        .fft_size = range_fft_size,
        .mean_removal_enabled = true,
        .window_config = {IFX_WINDOW_BLACKMANHARRIS, dev_config->num_samples_per_chirp, 0, 1},
        .is_normalized_window = 1};

    ifx_PPFFT_Config_t doppler_fft_config = (ifx_PPFFT_Config_t){
        .fft_type = IFX_FFT_TYPE_C2C,
        .fft_size = doppler_fft_size,
        .mean_removal_enabled = true,
        .window_config = {IFX_WINDOW_CHEBYSHEV, dev_config->num_chirps_per_frame, 100, 1},
        .is_normalized_window = 1};

    ifx_RDM_Config_t rdm_config = (ifx_RDM_Config_t){
        .spect_threshold = IFX_SPECT_THRESHOLD,
        .output_scale_type = IFX_SCALE_TYPE_LINEAR,
        .range_fft_config = range_fft_config,
        .doppler_fft_config = doppler_fft_config};

    ifx_DBF_Config_t dbf_config = (ifx_DBF_Config_t){
        .num_beams = 8,
        .num_antennas = 3,
        .min_angle = -45,
        .max_angle = 45,
        .d_by_lambda = 1000.0f};

    ifx_RAI_Config_t rai_config = (ifx_RAI_Config_t){
        .alpha_mti_filter = 0.5f,
        .dbf_config = dbf_config,
        .rdm_config = rdm_config,
        .num_antenna_array = 8,
        .num_of_images = 2};

    ifx_OSCFAR_Config_t oscfar_config = (ifx_OSCFAR_Config_t){
        .win_rank = 4,
        .guard_band = 2,
        .sample = 8,
        .pfa = 0.01,
        .coarse_scalar = 0.7};

    ifx_DBSCAN_Config_t dbscan_config = (ifx_DBSCAN_Config_t){
        .min_points = 2,
        .min_dist = 2.0f,
        .max_num_detections = 3,
    };

    app_context->rai_handle = ifx_rai_create(&rai_config);
    if ((ret = ifx_error_get()))
    {
        printf("app_radar_config: errorrr \r\n");
        return ret;
    }

    app_context->oscfar_handle = ifx_oscfar_create(&oscfar_config);
    app_context->dbscan_handle = ifx_dbscan_create(&dbscan_config);

    app_context->rai_cube_r = ifx_cube_create_r(dev_config->num_chirps_per_frame, dev_config->num_samples_per_chirp, dbf_config.num_beams);
    app_context->rdm_r = ifx_mat_create_r(range_fft_size / 2, doppler_fft_size);
    return (ifx_error_get());
}

ifx_Error_t app_radar_cleanup(app_handle_t *app_context)
{

    ifx_dbscan_destroy(app_context->dbscan_handle);
    ifx_oscfar_destroy(app_context->oscfar_handle);

    ifx_rai_destroy(app_context->rai_handle);
    ifx_cube_destroy_r(app_context->rai_cube_r);
    ifx_mat_destroy_r(app_context->rdm_r);
    return ifx_error_get();
}

ifx_Error_t app_radar_process(app_handle_t *app_context, ifx_Cube_R_t *frame)
{
    uint32_t col = IFX_CUBE_COLS(frame);
    uint32_t row = IFX_CUBE_ROWS(frame);
    uint32_t slice = IFX_CUBE_SLICES(frame);
    uint32_t rmax = 0;
    uint32_t smax = 0;
    uint32_t cmax = 0;

    static int abc = 0;
    abc++;

    printf("\r\n ");
    printf("frame data: col:%u  row:%u  slice:%u\r\n", col, row, slice);

    ifx_rai_run_r(app_context->rai_handle, frame, app_context->rai_cube_r);

    ifx_Float_t value = IFX_CUBE_AT(app_context->rai_cube_r, 0, 0, 0);
    ifx_Matrix_R_t rai_slice;

    for (uint32_t s = 0; s < IFX_CUBE_SLICES(app_context->rai_cube_r); s++)
    {
        uint32_t rm;
        uint32_t cm;

        ifx_cube_get_slice_r(app_context->rai_cube_r, s, &rai_slice);
        rdm_peak_search(&rai_slice, &rm, &cm);

        if (abc==1)
        {
            for (uint32_t col = 0; col < IFX_MAT_COLS(&rai_slice); col++)
            {
                ifx_Float_t value = IFX_MAT_AT(&rai_slice, rmax, col);
                printf("%10.6f  ", value);
            }
            printf("\r\n\n");
        }

        if (value < IFX_CUBE_AT(app_context->rai_cube_r, rm, cm, s))
        {
            value = IFX_CUBE_AT(app_context->rai_cube_r, rm, cm, s);
            rmax = rm;
            cmax = cm;
            smax = s;
        }
    }
    printf("rai_slice  rmax:%u  cmax:%u   slcie:%u    ..value:%10.6f\r\n", rmax, cmax, smax, value);

    ifx_Cube_C_t *rdm_cube_r = ifx_rai_get_range_doppler(app_context->rai_handle);

    uint32_t doppler_fft_size = IFX_CUBE_COLS(rdm_cube_r);
    uint32_t range_fft_size = IFX_CUBE_ROWS(rdm_cube_r);
    uint32_t num_virtual_antennas = IFX_CUBE_SLICES(rdm_cube_r);

    ifx_Matrix_C_t rdm_c;
    ifx_cube_get_slice_c(rdm_cube_r, 0, &rdm_c);
    for (uint32_t i = 0; i < IFX_MAT_ROWS(app_context->rdm_r); ++i)
    {
        ifx_Vector_C_t rdm_view;
        ifx_mat_get_rowview_c(&rdm_c, i, &rdm_view);

        ifx_Vector_R_t output_vec;
        ifx_mat_get_rowview_r(app_context->rdm_r, i, &output_vec);

        /* compute squared norm of spectrum */
        ifx_vec_abs_c(&rdm_view, &output_vec);

        // /* convert to linear or to dB */
        // spectrum2_to_linear(&output_vec, IFX_SPECT_THRESHOLD);
    }
    // do peak search

    rdm_peak_search(app_context->rdm_r, &rmax, &cmax);

    printf("range doppler doppler_fft_size=%u   range_fft_size=%u   num_virtual_antennas=%u   rmax:%u  cmax:%u\r\n",
           doppler_fft_size, range_fft_size, num_virtual_antennas, rmax, cmax);

    ifx_Vector_R_t *snr_vec_r = ifx_rai_get_snr(app_context->rai_handle);

    if (abc == 1)
    {
        printf("Getter function to access SNR result. signal-to-noise ratios \r\n");
        for (uint32_t ii = 0; ii < IFX_VEC_LEN(snr_vec_r); ii++)
        {
            printf("%10.6f  ", IFX_VEC_AT(snr_vec_r, ii));
        }
        printf("\r\n\n");
    }
    return (ifx_error_get());
}

int main(int argc, char **argv)
{
    app_t app_radar = {0};
    app_handle_t app_context = {0};
    int exitcode = 0;

    static const char *app_description = "App Radar";
    static const char *app_epilog = NULL;

    app_radar.app_description = app_description;
    app_radar.app_epilog = app_epilog;

    app_radar.app_init = (void *)&app_radar_init;
    app_radar.app_config = (void *)&app_radar_config;
    app_radar.app_process = (void *)&app_radar_process;
    app_radar.app_cleanup = (void *)&app_radar_cleanup;

    app_radar.default_metrics = &default_metrics;

    exitcode = app_start(argc, argv, &app_radar, &app_context);
    return exitcode;
}
