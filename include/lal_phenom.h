#ifndef LAL_PHENOM_H
#define LAL_PHENOM_H

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>

void generatePhenomLAL(
    const mass_t        mass_1, 
    const mass_t        mass_2, 
    const float64_t     sample_rate_hertz, 
    const timeUnit_t    duration, 
    const float64_t     inclination, 
    const length_t      distance, 
          float64_2_t **ret_strain
    ) {
	
	const int32_t num_samples = 
		(int32_t)floor(sample_rate_hertz*duration.seconds);
	
	REAL8TimeSeries *hplus  = NULL;
	REAL8TimeSeries *hcross = NULL;
	
	REAL8 S1x          = 0.0;
	REAL8 S1y          = 0.0;
	REAL8 S1z          = 0.0;
	REAL8 S2x          = 0.0;
	REAL8 S2y          = 0.0;
	REAL8 S2z          = 0.0;
	
	REAL8 phiRef       = 0.0;
	REAL8 longAscNodes = 0.0;
	REAL8 eccentricity = 0.0;
	REAL8 meanPerAno   = 0.0;
	REAL8 deltaT       = 1.0/sample_rate_hertz;
	REAL8 f_min        = 
		calcMinimumFrequency(
			mass_1, 
			mass_2, 
			duration
		);
	
	REAL8        f_ref       = 0.0;
	LALDict     *extraParams = NULL;
	Approximant  approximant = IMRPhenomXPHM;
	
	//Converting to SI:
		
	XLALSimInspiralTD(
		&hplus,
		&hcross,
		mass_1.kilograms,
		mass_2.kilograms,
		S1x,
		S1y,
		S1z,
		S2x,
		S2y,
		S2z,
		distance.meters,
		inclination,
		phiRef,
		longAscNodes,
		eccentricity,
		meanPerAno,
		deltaT,
		f_min,
		f_ref,
		extraParams,
		approximant
	);
	
	const int32_t waveform_num_samples = (int32_t)hplus->data->length;
	
	if (waveform_num_samples < num_samples) 
	{	
		fprintf(
			stderr, 
			"Warning! LAL Simulation not generating waveforms of desired "
			"num_samples.\n"
		);
	}
	
	size_t new_array_size = (size_t)num_samples * sizeof(float64_2_t);

	float64_2_t *strain = (float64_2_t*)malloc(new_array_size);
	int32_t new_waveform_index = 0;
    for (int32_t index = 0; index < num_samples; index++) 
	{	
		new_waveform_index = waveform_num_samples - num_samples - 1 + index;
		strain[index].x = (float64_t)hcross->data->data[new_waveform_index];
		strain[index].y = (float64_t)hplus->data->data[new_waveform_index];
    }
	
	free(hcross->data->data); free(hplus->data->data);

	*ret_strain = strain;
}

#endif