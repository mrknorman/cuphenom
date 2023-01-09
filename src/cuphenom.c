#include <math.h>
#include <inttypes.h>

#include "py_tools.h"

#define NUM_POLARIZATION_STATES 2

typedef struct {
	float x; 
	float y;
} float32_2_t;

#include "lal_phenom.h"


void addLinearArray(
          float   *array, 
    const float    min, 
    const float    max, 
    const int32_t  num_elements
    ) {
	
	// Calculate the increase in value between two subsequent array indices:
	const float step_increase = (max - min) / (float)num_elements; 

	for (int32_t index = 0; index < num_elements; index++) 
    {	
		array[index] = min + ((float)index)*step_increase;
	}
}

int32_t plotWaveform(
    const int32_t      verbosity,
	const float        duration_seconds,
          float32_2_t *strain,
	const int32_t      num_samples,
    const char        *output_file_name
    ) {
	
	const size_t array_size = sizeof(float)*(size_t)num_samples;
	
	float *h_cross = malloc(array_size);
	float *h_plus  = malloc(array_size);
	float *duration_array = malloc(array_size); 

	for (int32_t index = 0; index < num_samples; index++)
	{
		h_cross[index] = strain[index].x;
		h_plus[index] = strain[index].y;
	}	        
    addLinearArray(
        duration_array, 
        0.0f, 
        duration_seconds, 
        num_samples
    );

    // Setup axis:
    int32_t num_axis = NUM_POLARIZATION_STATES;
    axis_s y_axis = {
        .name = "strain",
        .label = "Strain"
    };
    axis_s x_axis = {
        .name = "time",
        .label = "Time (Seconds)"
    };
    axis_s axis_array[] = {
        x_axis,
        y_axis
    };

    // Setup axes:
	axes_s h_cross_axes = {
        .name     = "h_cross",
        .num_axis = num_axis,        
        .axis     = axis_array
    };
    axes_s h_plus_axes = {
        .name     = "h_plus",
        .num_axis = num_axis,        
        .axis     = axis_array
    };

    // Setup series:
    series_values_s duration_values = {
        .axis_name = "time",
        .values    = duration_array
    };
    series_values_s h_cross_values = {
        .axis_name = "strain",
        .values    = h_cross
    };
    series_values_s h_plus_values = {
        .axis_name = "strain",
        .values    = h_plus
    };
    
    series_values_s h_cross_series_values[] = {
        duration_values,
        h_cross_values
    };
	series_values_s h_plus_series_values[] = {
        duration_values,
        h_plus_values
    };

    series_s h_cross_series = {
        .label         = "h_cross",
        .axes_name     = "h_cross",
        .num_elements  = num_samples,
        .num_axis      = num_axis,
        .values        = h_cross_series_values
    };
    series_s h_plus_series = {
        .label         = "h_plus",
        .axes_name     = "h_plus",
        .num_elements  = num_samples,
        .num_axis      = num_axis,
        .values        = h_plus_series_values
    };

    axes_s axes[] = {
        h_cross_axes,
        h_plus_axes
    };
    series_s series[] = {
        h_cross_series,
        h_plus_series
    };

    figure_s figure = {
        .axes       = axes,
        .num_axes   = NUM_POLARIZATION_STATES,
        .series     = series,
        .num_series = NUM_POLARIZATION_STATES,
    };
	
    plotFigure(
        verbosity,
        figure,
        output_file_name
    ); 
    
    free(duration_array);
	free(h_plus);
	free(h_cross);
	
	return 0;
}

int32_t main(){

	const float   mass_1_msun       =   10.0f;
    const float   mass_2_msun       =   10.0f;
    const float   sample_rate_hertz = 8192.0f;
    const float   duration_seconds  =    1.0f;
    const float   inclination       =    0.0f;
    const float   distance_mpc      =   10.0f;
	
	float32_2_t *strain = NULL;
	const int32_t num_samples = 
		(int32_t)floorf(sample_rate_hertz*duration_seconds);
	
	generateLALInspiral(
		mass_1_msun, 
		mass_2_msun, 
		sample_rate_hertz, 
		num_samples, 
		inclination, 
		distance_mpc, 
		&strain
    );
	
	const char *output_file_name = 
		"../cuphenom_outputs/waveform_tests/lal_waveform";
		
	plotWaveform(
		STANDARD,
		duration_seconds,
        strain,
		num_samples,
		output_file_name
    );
		
	free(strain);
	return 0;
}