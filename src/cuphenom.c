#include <math.h>
#include <inttypes.h>
#include <time.h>

#include "py_tools.h"
#include "py_tools.h"
#include "units.h"

float64_t calcMinimumFrequency(
    const mass_t     mass_1,   //<-- Mass of first object.
    const mass_t     mass_2,   //<-- Mass of secondary object.
    const timeUnit_t duration //<-- Duration of signal.
) {
    /*
     * Calculates minimum frequency based on inputted masses.
     */
    
    const float64_t MC  = 
		pow(
			(((mass_1.kilograms*mass_2.kilograms)*(mass_1.kilograms*mass_2.kilograms)*
			  (mass_1.kilograms*mass_2.kilograms))/(mass_1.kilograms+mass_2.kilograms)),
			(1.0/5.0));
    float64_t fgw = (pow((duration.seconds/5.0),(-3.0/8.0)))*(1.0/(8.0*M_PI))
		*(pow((G_SI*MC/(C_SI*C_SI*C_SI)),(-5.0/8.0)));
  
    fgw = (1.0 > fgw) + (1.0 <= fgw)*fgw;
    
    return fgw;
}

#include "lal_phenom.h"
#include "cuda_phenom.h"

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
	const timeUnit_t   duration,
          float64_2_t *strain,
	const int32_t      num_samples,
    const char        *output_file_name
    ) {
	
	const size_t array_size = sizeof(float)*(size_t)num_samples;
	
	float *h_cross = malloc(array_size);
	float *h_plus  = malloc(array_size);
	float *duration_array = malloc(array_size); 

	for (int32_t index = 0; index < num_samples; index++)
	{
		h_cross[index] = (float)strain[index].x;
		h_plus[index]  = (float)strain[index].y;
	}	        
    addLinearArray(
        duration_array, 
        0.0f, 
        (float)duration.seconds, 
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

int32_t testRunTime(
	const mass_t     mass_1,
    const mass_t     mass_2,
    const float64_t  sample_rate_hertz,
    const timeUnit_t duration,
    const float64_t  inclination,
    const length_t   distance,
	const int32_t    num_tests
	) {
	
	float64_2_t *strain = NULL;
	
	clock_t start, end;
	float64_t execution_time_lal = 0.0, execution_time_cuda = 0.0;

	start = clock();
	for (int32_t index = 0; index < num_tests; index++)
	{
		generatePhenomLAL(
			mass_1, 
			mass_2, 
			sample_rate_hertz, 
			duration, 
			inclination, 
			distance, 
			&strain
		);
		free(strain);
	}
	
	end = clock();
	execution_time_lal = ((float64_t)(end - start))/CLOCKS_PER_SEC;
	
	start = clock();
	// CUDA:
	for (int32_t index = 0; index < num_tests; index++)
	{
		generatePhenomCUDA(
			mass_1,
			mass_2, 
			sample_rate_hertz, 
			duration, 
			inclination, 
			distance, 
			&strain
		);
			
		free(strain);
	}
	
	end = clock();
	execution_time_cuda = ((float64_t)(end - start))/CLOCKS_PER_SEC;
	
	printf("LAL: %f, CUDA: %f \n", execution_time_lal, execution_time_cuda);
	
	return 0;
}

int32_t main(){
	const float64_t   mass_1_msun       =   10.0f;
    const float64_t   mass_2_msun       =   10.0f;
    const float64_t   sample_rate_hertz = 8192.0f;
    const float64_t   duration_seconds  =    1.0f;
    const float64_t   inclination       =    0.0f;
    const float64_t   distance_mpc      = 1000.0f;
	
	float64_2_t *lal_strain = NULL, *cuda_strain = NULL;
		
	const mass_t mass_1 = initMassSolarMass(mass_1_msun);
	const mass_t mass_2 = initMassSolarMass(mass_2_msun);
	
	const length_t   distance = initLengthMpc(distance_mpc);
	const timeUnit_t duration = {.seconds = duration_seconds};
	
	const int32_t num_samples = 
		(int32_t)floor(sample_rate_hertz*duration.seconds);
	
	// LAL:
	generatePhenomLAL(
		mass_1, 
		mass_2, 
		sample_rate_hertz, 
		duration, 
		inclination, 
		distance, 
		&lal_strain
	);
	
	generatePhenomCUDA(
		mass_1,
		mass_2, 
		sample_rate_hertz, 
		duration, 
		inclination, 
		distance, 
		&cuda_strain
    );
	
	float64_2_t *difference = malloc(sizeof(float64_2_t)*(size_t)num_samples);
	float64_2_t sum = {.x = 0.0f, .y = 0.0f};
	for (int32_t index = 0; index < num_samples; index++)
	{
		difference[index].x = fabs(lal_strain[index].x - cuda_strain[index].x);
		difference[index].y = fabs(lal_strain[index].y - cuda_strain[index].y);
		
		sum.x += difference[index].x; sum.y += difference[index].y;
	}
	
	printf("SUM: %.4e, %.4e \n", sum.x, sum.y);
	
	// Plotting:
		
	char *output_file_name = 
		"../cuphenom_outputs/waveform_tests/lal_waveform";
		
	plotWaveform(
		STANDARD,
		duration,
        lal_strain,
		num_samples,
		output_file_name
    );
	
	output_file_name =  "../cuphenom_outputs/waveform_tests/cuda_waveform";
	
	plotWaveform(
		STANDARD,
		duration,
        cuda_strain,
		num_samples,
		output_file_name
    );
	
	free(lal_strain); free(cuda_strain);
	
	testRunTime(
		mass_1,
		mass_2,
		sample_rate_hertz,
		duration,
		inclination,
		distance,
		1000
	);
	
	return 0;
}