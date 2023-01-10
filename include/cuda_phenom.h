#ifndef CUDA_PHENOM_H
#define CUDA_PHENOM_H

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>


typedef struct {
	mass_t    _1;
	mass_t    _2;
	mass_t    total;
	mass_t    reduced;
	float64_t symmetric_ratio;
} mass_properties_t;

mass_properties_t initMassProperties(
	const mass_t mass_1,
	const mass_t mass_2
	) {
	
	mass_properties_t mass;
	
	mass._1              = mass_1;
	mass._2              = mass_2;
	mass.total           = addMasses(mass._1, mass._2);
    mass.reduced = 
		divideMasses(multiplyMasses(mass._1, mass._2), mass.total);
    mass.symmetric_ratio = mass.reduced.kilograms / mass.total.kilograms;
	
	return mass;
}

typedef struct {
	float64_t x;
	float64_t y;
	float64_t z;
} spin_t;

float64_t calculateSpinNorm(
	const spin_t spin
	) {
	
	return spin.x*spin.x + spin.y*spin.y + spin.z*spin.z;
}

typedef struct {
	spin_t    _1;
	spin_t    _2;
} spin_properties_t;

spin_properties_t initSpinProperties(
	const float64_t *spin_1,
	const float64_t *spin_2
	) {
	
	spin_properties_t spin;
	
	spin._1 = (spin_t) {
		.x = spin_1[0],
		.y = spin_1[1],
		.z = spin_1[2]
	};
	
	spin._2 = (spin_t) {
		.x = spin_2[0],
		.y = spin_2[1],
		.z = spin_2[2]
	};

	return spin;
}

typedef struct {
	mass_properties_t mass;
	spin_properties_t spin;
	length_t          distance;
} binarySystemProperties_t;

inline float64_t InspiralTaylorT2Timing_0PNCoeff(
	const mass_t total_mass,
    const float64_t sym_mass_ratio
) {
    return -5.0*total_mass.seconds/(256.0*sym_mass_ratio);
}

inline float64_t InspiralTaylorT2Timing_2PNCoeff(
	 const float64_t sym_mass_ratio)
{
	return 7.43/2.52 + 11./3. * sym_mass_ratio;
}

inline float64_t InspiralTaylorT2Timing_4PNCoeff(
	const float64_t sym_mass_ratio)
{
    return 30.58673/5.08032 + 54.29/5.04*sym_mass_ratio 
		 + 61.7/7.2*sym_mass_ratio*sym_mass_ratio;
}
  
timeUnit_t InspiralChirpTimeBound(
	const float64_t         fstart, 
	const mass_properties_t mass,
	const spin_properties_t spin
) {
	 // over-estimate of chi
	 const float64_t chi = fabs(
		    ((fabs(spin._1.z) > fabs(spin._2.z)) * spin._1.z)
		 +  ((fabs(spin._1.z) <= fabs(spin._2.z)) * spin._2.z)
		 );
	 
     const float64_t c0 = 
		 fabs(InspiralTaylorT2Timing_0PNCoeff(mass.total, mass.symmetric_ratio));
     const float64_t c2 = 
		 InspiralTaylorT2Timing_2PNCoeff(mass.symmetric_ratio);
	
     /* the 1.5pN spin term is in TaylorT2 is 8*beta/5
      * where beta = (113/12 + (25/4)(m2/m1))*(s1*m1^2/M^2) + 2 <-> 1
      * [Cutler & Flanagan, Physical Review D 49, 2658 (1994), Eq. (3.21)]
      * which can be written as (113/12)*chi - (19/6)(s1 + s2)
      * and we drop the negative contribution 
	  */
     const float64_t c3 = (226.0/15.0) * chi;
	 
     /* there is also a 1.5PN term with eta, but it is negative so do not 
	  * include it */
     const float64_t c4 = InspiralTaylorT2Timing_4PNCoeff(mass.symmetric_ratio);
     const float64_t v = cbrt(M_PI * G_SI * mass.total.kilograms * fstart) / C_SI;
	 
     return initTimeSeconds(c0 * pow(v, -8) * (1.0 + (c2 + (c3 + c4 * v) * v) * v * v));
}

float64_t InspiralFinalBlackHoleSpinBound(
	const spin_properties_t spin
	) {
     const float64_t maximum_black_hole_spin = 0.998;
     /* lower bound on the final plunge, merger, and ringdown time here the
      * final black hole spin is overestimated by using the formula in Tichy and
      * Marronetti, Physical Review D 78 081501 (2008), Eq. (1) and Table 1, for
      * equal mass black holes, or the larger of the two spins (which covers the
      * extreme mass case) 
	  */
	  
     float64_t final_spin_upper_bound = 0.686 + 0.15 * (spin._1.z + spin._2.z);
	 final_spin_upper_bound = 
		   (final_spin_upper_bound < fabs(spin._1.z))*fabs(spin._1.z) 
		 + (final_spin_upper_bound > fabs(spin._1.z))*final_spin_upper_bound;
	 final_spin_upper_bound = 
		   (final_spin_upper_bound < fabs(spin._2.z))*fabs(spin._2.z) 
		 + (final_spin_upper_bound > fabs(spin._2.z))*final_spin_upper_bound;

     /* it is possible that |S1z| or |S2z| >= 1, but s must be less than 1
     * (0th law of thermodynamics) so provide a maximum value for s */	
	 final_spin_upper_bound = 
			(final_spin_upper_bound > maximum_black_hole_spin)
				* maximum_black_hole_spin
		  + (final_spin_upper_bound < maximum_black_hole_spin)
			    * final_spin_upper_bound;

     return final_spin_upper_bound;
}

inline timeUnit_t InspiralMergeTimeBound(
	const mass_properties_t mass
) {
	return initTimeSeconds(2.0*M_PI*((9.0*mass.total.meters)/(C_SI/3.0)));
}

timeUnit_t InspiralRingdownTimeBound(
	const mass_properties_t mass,
	const float64_t         final_spin_upper_bound
	) {
     const float64_t nefolds = 11; /* Waveform generators only go up to 10 */
  
     /* These values come from Table VIII of Berti, Cardoso, and Will with n=0, m=2 */
	 const float64_t f[] = {1.5251, -1.1568,  0.1292}; 
	 const float64_t q[] = {0.7000,  1.4187, -0.4990}; 

     const float64_t omega = (f[0] + f[1] * pow(1.0 - final_spin_upper_bound, f[2])) 
	                   / mass.total.seconds;
     const float64_t Q = q[0] + q[1] * pow(1.0 - final_spin_upper_bound, q[2]);
     const float64_t tau = 2.0 * Q / omega; /* see Eq. (2.1) of Berti, Cardoso, and Will */
	 
	 return initTimeSeconds(nefolds * tau);
}

inline double InspiralTaylorT3Frequency_0PNCoeff(
         const mass_t mass
	) {	
    return 1.0 / (8.0 * M_PI * mass.seconds);
}

float64_t InspiralChirpStartFrequencyBound(
	const timeUnit_t        chirp_time, 
	const mass_properties_t mass
	) {
	 
	double c0 = InspiralTaylorT3Frequency_0PNCoeff(mass.total);
	return c0 * pow(5.0 * mass.total.seconds / (mass.symmetric_ratio * chirp_time.seconds), 3.0 / 8.0);
}

mass_properties_t initRedshiftedMass(
	      mass_t    mass_1,
	      mass_t    mass_2,
	const float64_t redshift
	) {
	
	 // Apply redshift correction to dimensionful source-frame quantities:
	 mass_1   = scaleMass(mass_1, 1.0 + redshift);
	 mass_2   = scaleMass(mass_2, 1.0 + redshift);
	  
	// Create mass_properties struct to hold common mass properties:
	return initMassProperties(mass_1, mass_2);
}

//#include <choose.h>

static void checkSystemParameters(
	const mass_properties_t mass,
	const spin_properties_t spin, 
	const timeUnit_t        deltaT,
	const float64_t         f_min
	) {

	if (deltaT.seconds > 1.0)
	{
		fprintf(
			stderr,
			"Warning - %s: Large value of deltaT = %e (s). requested.\nPerhaps sample"
			" rate and time step size were swapped?\n", 
			__func__, deltaT.seconds
		);
	}
	if (deltaT.seconds < 1.0/16385.0)
	{
		fprintf(
			stderr,
			"Warning - %s: Small value of deltaT = %e (s). requested.\nCheck "
			" for errors, this could create very large time series.\n", 
			__func__, deltaT.seconds
		);
	}
	if (mass._1.kilograms < initMassSolarMass(0.09).kilograms)
	{
		fprintf(
			stderr, 
			"Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\n"
			" Perhaps you have a unit conversion error?\n", 
			__func__, mass._1.kilograms, mass._1.kilograms/ MASS_SUN_KILOGRAMS
		);
	}
	if (mass._2.kilograms < initMassSolarMass(0.09).kilograms)
	{
         fprintf(stderr,
			"Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\n"
			"Perhaps you have a unit conversion error?\n", 
			__func__, mass._2.kilograms, mass._2.kilograms/MASS_SUN_KILOGRAMS);
	}
    if (mass.total.kilograms > initMassSolarMass(1000).kilograms)
	{
         fprintf(
			 stderr,
			"Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun)"
			" requested.\nSignal not likely to be in band of ground-based"
			" detectors.\n", 
			__func__, 
			mass.total.kilograms, 
			mass.total.kilograms/MASS_SUN_KILOGRAMS
		);
	}
    if (calculateSpinNorm(spin._1) > 1.000001)
	{
         fprintf(
			 stderr, 
			 "Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested.\nAre you "
			 "sure you want to violate the Kerr bound?\n", 
			 __func__, spin._1.x, spin._1.y, spin._1.z
		);
	}
    if(calculateSpinNorm(spin._2) > 1.000001)
	{
		fprintf(
			stderr, 
			"Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested.\nAre you "
			"sure you want to violate the Kerr bound?\n", 
			__func__, spin._2.x, spin._2.y, spin._2.z
		);
	}
    if(f_min < 1.0)
	{
		fprintf(
			stderr, 
			"Warning - %s: Small value of fmin = %e requested.\nCheck for "
			"errors, this could create a very long waveform.\n", 
			__func__, f_min);
	}
    if(f_min > 40.000001)
	{
         fprintf(
			 stderr, 
			 "Warning - %s: Large value of fmin = %e requested.\nCheck for "
			 "errors, the signal will start in band.\n", 
			 __func__, f_min
		);
	}
}

int32_t InspiralTDFromTD(
     REAL8TimeSeries **hplus,  // <-- +-polarization waveform */
     REAL8TimeSeries **hcross, // <-- x-polarization waveform */
     mass_t mass_1,            // <-- mass of companion 1 (mass_t) */
     mass_t mass_2,            // <-- mass of companion 2 (mass_t) */
	 float64_t *spin_1,        // <-- vector of dimensionless spins of companion 1 {x,y,z}
	 float64_t *spin_2,        // <-- vector of dimensionless spins of companion 2 {x,y,z}
     length_t distance,        // <-- distance of source (m) */
     REAL8 inclination,        // <-- inclination of source (rad) */
     REAL8 phiRef,             // <-- reference orbital phase (rad) */
     REAL8 longAscNodes,       // <-- longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
     REAL8 eccentricity,       // <-- eccentrocity at reference epoch */
     REAL8 meanPerAno,         // <-- mean anomaly of periastron */
     timeUnit_t deltaT,        // <-- sampling interval (timeUnit_t) */
     REAL8 f_min,              // <-- starting GW frequency (Hz) */
     REAL8 f_ref,              // <-- reference GW frequency (Hz) */
	 REAL8 redshift,           // r--edshift correction to dimensionful source-frame quantities 
     LALDict *LALparams,       // <-- LAL dictionary containing accessory parameters */
     Approximant approximant   // <-- post-Newtonian approximant to use for waveform production */
) {
	
	// Hard coded constants:
	const float64_t extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
    const float64_t extra_cycles        = 3.0; /* more extra time measured in cycles at the starting frequency */
	
	const spin_properties_t spin = initSpinProperties(spin_1, spin_2);
	const mass_properties_t mass = initRedshiftedMass(mass_1, mass_2, redshift);
	distance = scaleLength(distance, 1.0 + redshift);  /* change from comoving (transverse) distance to luminosity distance */
	
	const float64_t original_f_min = f_min; /* f_min might be overwritten below, so keep original value */
	
	// Because we are using precessing waveform set f_ref to f_min:
	f_ref = f_min; 

	// If the requested low frequency is below the lowest Kerr ISCO
	// frequency then change it to that frequency:
	float64_t kerr_isco_frequency = 1.0 / (pow(9.0, 1.5) * M_PI * mass.total.seconds);
	if (f_min > kerr_isco_frequency)
		f_min = kerr_isco_frequency;
	
	// General sanity check the input parameters - only give warnings! 
	checkSystemParameters(
		mass,
		spin, 
		deltaT,
		f_min
	);
	
	// Upper bound on the chirp time starting at f_min:
	const timeUnit_t chirp_time_upper_bound = 
		InspiralChirpTimeBound(f_min, mass, spin);
	
	// Upper bound on the final black hole spin:
     const float64_t final_spin_upper_bound = 
		 InspiralFinalBlackHoleSpinBound(spin);
	
    // Upper bound on the final plunge, merger, and ringdown time:
	const timeUnit_t plunge_time_upper_bound = 
		addTimes(
			2,
			InspiralMergeTimeBound(mass),
			InspiralRingdownTimeBound(mass, final_spin_upper_bound)
		);
  
	// Extra time to include for all waveforms to take care of situations
	// where the frequency is close to merger (and is sweeping rapidly):
	// this is a few cycles at the low frequency:
	const timeUnit_t extra_time = initTimeSeconds(extra_cycles / f_min);
	 
	// time domain approximant: condition by generating a waveform
	// with a lower starting frequency and apply tapers in the
	// region between that lower frequency and the requested
	// frequency f_min; here compute a new lower frequency:
	 
	timeUnit_t region_max_duration = 		 
		addTimes(
			 3,
			 chirp_time_upper_bound,
			 plunge_time_upper_bound,
			 extra_time
		 );
	
	region_max_duration = scaleTime(
		region_max_duration,
		(1.0 + extra_time_fraction)
	);
	
	float64_t fstart = 
		 InspiralChirpStartFrequencyBound(
			 region_max_duration, 
			 mass
		 );
  
     // Generate the waveform in the time domain starting at fstart :
     int32_t retval = 
		 XLALSimInspiralChooseTDWaveform(
			 hplus, 
			 hcross, 
			 mass._1.kilograms, 
			 mass._2.kilograms, 
			 spin._1.x, spin._1.y, spin._1.z, 
			 spin._2.x, spin._2.y, spin._2.z, 
			 distance.meters, 
			 inclination, 
			 phiRef, 
			 longAscNodes, 
			 eccentricity, 
			 meanPerAno, 
			 deltaT.seconds,
			 fstart, 
			 f_ref, 
			 LALparams, 
			 approximant
	);
	
     if (retval < 0)
         XLAL_ERROR(XLAL_EFUNC);
  
     /* condition the time domain waveform by tapering in the extra time
         * at the beginning and high-pass filtering above original f_min */
     XLALSimInspiralTDConditionStage1(
		 *hplus, 
		 *hcross, 
		 extra_time_fraction * chirp_time_upper_bound.seconds + extra_time.seconds, 
		 original_f_min
	 );
  
     return 0;
 }


void generatePhenomCUDA(
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
	
	float64_t spin_1[] = {0.0, 0.0, 0.0};  
	float64_t spin_2[] = {0.0, 0.0, 0.0};  
	
	REAL8 phiRef       = 0.0;
	REAL8 longAscNodes = 0.0;
	REAL8 eccentricity = 0.0;
	REAL8 meanPerAno   = 0.0;
	timeUnit_t deltaT  = initTimeSeconds(1.0/sample_rate_hertz);
	REAL8 f_min        = 
		calcMinimumFrequency(
			mass_1, 
			mass_2, 
			duration
		);
	
	REAL8 f_ref       = 0.0;
	REAL8 redshift    = 0.0;

	LALDict     *extraParams = NULL;
	Approximant  approximant = IMRPhenomXPHM;
	
	InspiralTDFromTD(
		&hplus,
		&hcross,
		mass_1,
		mass_2,
		spin_1,
		spin_2,
		distance,
		inclination,
		phiRef,
		longAscNodes,
		eccentricity,
		meanPerAno,
		deltaT,
		f_min,
		f_ref,
		redshift, 
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
		strain[index].x = (float)hcross->data->data[new_waveform_index];
		strain[index].y = (float)hplus->data->data[new_waveform_index];
    }
	
	free(hcross->data->data); free(hplus->data->data);

	*ret_strain = strain;
}

#endif