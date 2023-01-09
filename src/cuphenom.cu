#include <inttypes.h>
#include "lal_phenom.h"

int32_t main(){

	const float   mass_1_msun       =   10.0f;
    const float   mass_2_msun       =   10.0f;
    const float   sample_rate_hertz = 8192.0f;
    const int32_t num_samples       = 8192.0f;
    const float   inclination       =    0.0f;
    const float   distance_mpc      =   10.0f;
	
	float2 *strain = NULL;
	
	generateLALInspiral(
		mass_1_msun, 
		mass_2_msun, 
		sample_rate_hertz, 
		num_samples, 
		inclination, 
		distance_mpc, 
		&strain
    );
	
	return 0;
}