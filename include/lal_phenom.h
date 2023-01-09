#ifndef DATASET_WAVEFORM_H
#define DATASET_WAVEFORM_H

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%   Hardwired constants.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

const float G           = 6.673f*1e-11f;       //-- the gravitational constant in m^3/(kg*s^2)
const float HPLANK      = 6.626068e-34f;      //-- Planck's constant in m^2*kg/s
const float C           = 299792458.0f;         //-- the speed of light in m/s
const float PC_METERS   = 3.0856775807f*1e16f; //-- a parsec in m
const float MASS_SUN_KG = 1.98892f*1e30f;      //-- solar mass in kg

float calcMinimumFrequency(
    const float mass_1_msun,   //<-- Mass of first object.
    const float mass_2_msun,   //<-- Mass of secondary object.
    const float duration_seconds  //<-- duration_seconds of signal.
) {
    /*
     * Calculates minimum frequency based on inputted masses.
     */
    
    const float MC  = 
		powf(
			(((mass_1_msun*mass_2_msun)*(mass_1_msun*mass_2_msun)*(mass_1_msun*mass_2_msun))/(mass_1_msun+mass_2_msun)),
			(1.0f/5.0f))*MASS_SUN_KG;
    float fgw = (powf((duration_seconds/5.0f),(-3.0f/8.0f)))*(1.0f/(8.0f*(float)M_PI))
		*(powf((G*MC/(C*C*C)),(-5.0f/8.0f)));
  
    fgw = (1.0 > fgw) + (1.0 <= fgw)*fgw;
    
    return fgw;
}

void generateLALInspiral(
    const float         mass_1_msun, 
    const float         mass_2_msun, 
    const float         sample_rate_hertz, 
    const int32_t       num_samples, 
    const float         inclination, 
    const float         distance_mpc, 
          float32_2_t **ret_strain
    ) {
    
    const float duration_seconds = (float)num_samples/sample_rate_hertz;
	 
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
			mass_1_msun, 
			mass_2_msun, 
			duration_seconds
		);
	
	REAL8        f_ref       = 0.0;
	LALDict     *extraParams = NULL;
	Approximant  approximant = IMRPhenomD;
	
	//Converting to SI:
	
	REAL8 mass_1_kg = mass_1_msun*MASS_SUN_KG;
    REAL8 mass_2_kg = mass_2_msun*MASS_SUN_KG;
	
	REAL8 distance_meters = distance_mpc*PC_METERS*10E6f;
	
	XLALSimInspiralTD(
		&hplus,
		&hcross,
		mass_1_kg,
		mass_2_kg,
		S1x,
		S1y,
		S1z,
		S2x,
		S2y,
		S2z,
		distance_meters,
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
	
	size_t new_array_size = (size_t)num_samples * sizeof(float32_2_t);

	float32_2_t *strain = (float32_2_t*)malloc(new_array_size);
	int32_t new_waveform_index = 0;
    for (int32_t index = 0; index < num_samples; index++) 
	{	
		new_waveform_index = waveform_num_samples - num_samples - 1 + index;
		strain[index].x = (float)hcross->data->data[new_waveform_index];
		strain[index].y = (float)hplus->data->data[new_waveform_index];
    }
	
	free(hcross->data->data); free(hplus->data->data);

	*ret_strain = strain;
}

#endif