#ifndef CUDA_PHENOM_H
#define CUDA_PHENOM_H

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include <cuda_maths.h>

//#include "phenomd.h"

int32_t InspiralTDFromTD(
     REAL8TimeSeries **hplus,  // <-- +-polarization waveform 
     REAL8TimeSeries **hcross, // <-- x-polarization waveform 
     mass_t mass_1,            // <-- mass of companion 1 (mass_t) 
     mass_t mass_2,            // <-- mass of companion 2 (mass_t) 
	 float64_t *spin_1,        // <-- vector of dimensionless spins of companion 1 {x,y,z}
	 float64_t *spin_2,        // <-- vector of dimensionless spins of companion 2 {x,y,z}
     length_t distance,        // <-- distance of source (m) 
     REAL8 inclination,        // <-- inclination of source (rad) 
     REAL8 phiRef,             // <-- reference orbital phase (rad) 
     REAL8 longAscNodes,       // <-- longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation 
     REAL8 eccentricity,       // <-- eccentrocity at reference epoch 
     REAL8 meanPerAno,         // <-- mean anomaly of periastron 
     timeUnit_t deltaT,        // <-- sampling interval (timeUnit_t) 
     REAL8 f_min,              // <-- starting GW frequency (Hz) 
     REAL8 f_ref,              // <-- reference GW frequency (Hz) 
	 REAL8 redshift,           // r--edshift correction to dimensionful source-frame quantities 
     LALDict *LALparams,       // <-- LAL dictionary containing accessory parameters 
     Approximant approximant   // <-- post-Newtonian approximant to use for waveform production 
) {
	
	// Hard coded constants:
	
	// Fraction of waveform duration to add as extra time for tapering:
	const float64_t extra_time_fraction = 0.1;
	// More extra time measured in cycles at the starting frequency:
    const float64_t extra_cycles        = 3.0; 
	
	// Init property structures:
	const spin_properties_t spin = initSpinProperties(spin_1, spin_2);
	const mass_properties_t mass = initRedshiftedMass(mass_1, mass_2, redshift);
	
	// Change from comoving (transverse) distance to luminosity distance:
	distance = scaleLength(distance, 1.0 + redshift);  
	
	// f_min is overwritten below, so keep original value:
	const float64_t original_f_min = f_min; 
	
	// Because we are using precessing waveform set f_ref to f_min:
    
   	if (approximant == IMRPhenomXPHM) 
	{
        f_ref = f_min;
    } 

	// If the requested low frequency is below the lowest Kerr ISCO
	// frequency then change it to that frequency:
	float64_t kerr_isco_frequency = 
		1.0 / (pow(9.0, 1.5) * M_PI * mass.total.seconds);
	if (f_min > kerr_isco_frequency)
		f_min = kerr_isco_frequency;
	
	// Calculate f_max:
    double f_max = 0.5 / deltaT.seconds;

	// General sanity check the input parameters. This will only give warnings:
	checkSystemParameters(
		mass,
		spin, 
		deltaT,
		f_min
	);
		
	// Calculate time boundaries:

	// Upper bound on the chirp time starting at f_min with extra time to 
	// include for all waveforms to take care of situations where the frequency 
	// is close to merger (and is sweeping rapidly) this is a few cycles at the 
	// low frequency:
	const timeUnit_t extra_time = initTimeSeconds(extra_cycles / f_min);
	const timeUnit_t chirp_time_upper_bound = 
			addTimes(2, InspiralChirpTimeBound(f_min, mass, spin), extra_time);

	// Upper bound on the plunge and merger time:
	const timeUnit_t merge_time_upper_bound = InspiralMergeTimeBound(mass);

	// Upper bound on the ringdown time:
	const timeUnit_t ringdown_time_upper_bound = 
		InspiralRingdownTimeBound(mass, spin);
	
	const timeUnit_t total_time_upper_bound = 
		scaleTime(
			addTimes(
				3,
				chirp_time_upper_bound,
				merge_time_upper_bound,
				ringdown_time_upper_bound
			),
			(1.0 + extra_time_fraction)
		);
	
	// time domain approximant: condition by generating a waveform
	// with a lower starting frequency and apply tapers in the
	// region between that lower frequency and the requested
	// frequency f_min; here compute a new lower frequency:
	
	float64_t fstart = 
		InspiralChirpStartFrequencyBound(
			 total_time_upper_bound, 
			 mass
		);
	
	// generate the conditioned waveform in the frequency domain
    //note: redshift factor has already been applied above
    //set deltaF = 0 to get a small enough resolution
	
	// GPS times:
    timeUnit_t hplus_gps  = initTimeSeconds(0.0);
	timeUnit_t hcross_gps = initTimeSeconds(0.0);
	
	COMPLEX16FrequencySeries *hptilde = NULL;
	COMPLEX16FrequencySeries *hctilde = NULL;
	
	// total expected chirp length includes merger 
	size_t chirplen = (size_t) 
		round(total_time_upper_bound.seconds / deltaT.seconds);
	
	int32_t retval = 0;
	{
		float64_t deltaF = 0;
		float64_t new_fmin = fstart;
		
		// Apply condition that f_max rounds to the next power-of-two multiple
		// of deltaF.
		// Round f_max / deltaF to next power of two.
		// Set f_max to the new Nyquist frequency.
		// The length of the chirp signal is then 2 * f_nyquist / deltaF.
		// The time spacing is 1 / (2 * f_nyquist) 
		float64_t f_nyquist    = f_max;    
		int32_t   chirplen_exp =  0;

		if (deltaF != 0) 
		{
			int32_t n = (int32_t) round(f_max / deltaF);
			if ((n & (n - 1))) { // not a power of 2
				frexp(n, &chirplen_exp);
				f_nyquist = ldexp(1.0, chirplen_exp) * deltaF;
				fprintf(
					stderr,
					"f_max/deltaF = %g/%g = %g is not a power of two: changing "
					" f_max to %f", 
					f_max, deltaF, f_max/deltaF, f_nyquist
				);
			 }
		}
		double new_deltaT = 0.5 / f_nyquist;

		double new_tchirp = 
			InspiralChirpTimeBound(new_fmin, mass, spin).seconds;
		double new_tmerge = 
			merge_time_upper_bound.seconds + ringdown_time_upper_bound.seconds;
		
		double new_fstart = 
			XLALSimInspiralChirpStartFrequencyBound(
				(1.0 + extra_time_fraction) * new_tchirp, 
				mass._1.kilograms, 
				mass._2.kilograms
				);
		new_tchirp = InspiralChirpTimeBound(new_fstart, mass, spin).seconds;

		// need a long enough segment to hold a whole chirp with some padding
		// length of the chirp in samples
		chirplen = 
			(size_t)round((new_tchirp + new_tmerge + 2.0 * extra_time.seconds) 
			/ new_deltaT);
		
		// make chirplen next power of two 
		frexp((double)chirplen, &chirplen_exp);
		chirplen = (size_t)ldexp(1.0, chirplen_exp);

		if (deltaF == 0.0)
		{
			deltaF = 1.0 / ((double)chirplen * new_deltaT);
		}
		else if (deltaF > 1.0 / ((double)chirplen * new_deltaT))
		{
			fprintf(
				stderr, 
				" Specified frequency interval of %f Hz is too large for a"
				" chirp of duration %f s", 
				deltaF, 
				(double)chirplen * new_deltaT
			);
		}
		/*
		retval = SimInspiralFD(
				&hptilde, 
				&hctilde, 
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
				deltaF, 
				new_fstart, 
				f_max, 
				f_ref, 
				LALparams,
			    approximant
		);
		
		*/
		if (approximant == IMRPhenomXPHM) 
		{
			retval = XLALSimIMRPhenomXPHM(
				&hptilde, 
				&hctilde, 
				mass._1.kilograms, 
				mass._2.kilograms,
				spin._1.x, spin._1.y, spin._1.z, 
				spin._2.x, spin._2.y, spin._2.z, 
				distance.meters, 
				inclination,
				phiRef, 
				new_fstart, 
				f_max, 
				deltaF, 
				f_ref, 
				LALparams
			);
		}
		else if (approximant == IMRPhenomD)
		{
			double cfac = cos(inclination);
			double pfac = 0.5 * (1. + cfac*cfac);
			
			// Phenom D:
			retval = XLALSimIMRPhenomDGenerateFD(
				&hptilde, 
				phiRef, 
				f_ref, 
				deltaF, 
				mass._1.kilograms, 
				mass._2.kilograms,
				spin._1.z, 
				spin._2.z, 
				new_fstart, 
				f_max, 
				distance.meters, 
				LALparams, 
				4
			);
				
			//Produce both polarizations
			hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
				&(hptilde->epoch), hptilde->f0, hptilde->deltaF,
				&(hptilde->sampleUnits), hptilde->data->length);

			for(int32_t j = 0; j < hptilde->data->length; j++) 
			{
				hctilde->data->data[j] = -I*cfac * hptilde->data->data[j];
				hptilde->data->data[j] *= pfac;
			}
		}
		
		
		// revise (over-)estimate of chirp from new start frequency

		 // taper frequencies between fstart and f_min
		 int32_t k0 = round(new_fstart / hptilde->deltaF);
		 int32_t k1 = round(new_fmin / hptilde->deltaF);
		 
		 // make sure it is zero below fstart
		 for (int32_t k = 0; k < k0; ++k) 
		 {
			 hptilde->data->data[k] = 0.0;
			 hctilde->data->data[k] = 0.0;
		 }
		 // taper between fstart and f_min
		 for (int32_t k = k0; k < k1; ++k) 
		 {	 
			 double w = 0.5 - 0.5 * cos(M_PI * (k - k0) / (double)(k1 - k0));
			 hptilde->data->data[k] *= w;
			 hctilde->data->data[k] *= w;
		}
		// make sure Nyquist frequency is zero
		hptilde->data->data[hptilde->data->length - 1] = 0.0;
		hctilde->data->data[hctilde->data->length - 1] = 0.0;

		// we want to make sure that this waveform will give something
		// sensible if it is later transformed into the time domain:
		// to avoid the end of the waveform wrapping around to the beginning,
		// we shift waveform backwards in time and compensate for this
		// shift by adjusting the epoch 
		
		// Integer number of time samples:
		float64_t tshift = round(new_tmerge / new_deltaT) * new_deltaT; 

		 for (int32_t k = 0; k < hptilde->data->length; ++k) {
			 double complex phasefac = cexp(2.0 * M_PI * I * k * deltaF * tshift);
			 hptilde->data->data[k] *= phasefac;
			 hctilde->data->data[k] *= phasefac;
		 }
	}
		
	// we want to make sure that this waveform will give something
    // sensible if it is later transformed into the time domain:
    // to avoid the end of the waveform wrapping around to the beginning,
    // we shift waveform backwards in time and compensate for this
    // shift by adjusting the epoch -- note that XLALSimInspiralFD
    // uarantees that there is extra padding to do this 
	
	double deltaF = 1.0 / (double) ((int32_t) chirplen * deltaT.seconds);
	 
	// Integer number of samples 
	float64_t tshift = round(extra_time.seconds / deltaT.seconds) * deltaT.seconds; 
    for (size_t k = 0; k < hptilde->data->length; ++k) 
	{
		double complex phasefac = cexp(
			2.0 * M_PI * I * k * hptilde->deltaF * tshift);
		hptilde->data->data[k] *= phasefac;
		hctilde->data->data[k] *= phasefac;
    }
	hplus_gps  = addTimes(2, hplus_gps, deltaT);
	hcross_gps = addTimes(2, hcross_gps, deltaT);
	
	// Transform the waveform into the time domain:
	
	const LALUnit lalStrainUnit = 
		{ 0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
	
	size_t num_waveform_samples = 2 * (hptilde->data->length - 1);
	
	*hplus = XLALCreateREAL8TimeSeries(
		"H_PLUS", 
		&hptilde->epoch, 
		0.0, 
		deltaT.seconds,
		&lalStrainUnit, 
		num_waveform_samples
	);
	
	*hcross = XLALCreateREAL8TimeSeries(
		"H_CROSS", 
		&hctilde->epoch, 
		0.0, 
		deltaT.seconds, 
		&lalStrainUnit, 
		num_waveform_samples
	);
		
	performIRFFT64(
		hptilde->data->data,
		hctilde->data->data,
		(*hplus)->data->data,
		(*hcross)->data->data,
		deltaT,
		(int32_t)num_waveform_samples
	);
	
	// apply time domain filter at fstart
	XLALHighPassREAL8TimeSeries(*hplus, fstart, 0.99, 8);
	XLALHighPassREAL8TimeSeries(*hcross, fstart, 0.99, 8);

	// amount to snip off at the end is tshift 
	size_t end = 
		(*hplus)->data->length - (size_t) round(tshift / deltaT.seconds);

	// snip off extra time at beginning and at the end 
	int32_t snip_start = (int32_t) (end - chirplen);
	
	XLALResizeREAL8TimeSeries(*hplus, snip_start, chirplen);
	XLALResizeREAL8TimeSeries(*hcross, snip_start, chirplen);

	// clean up
	XLALDestroyCOMPLEX16FrequencySeries(hptilde);
	XLALDestroyCOMPLEX16FrequencySeries(hctilde);
	
    if (retval < 0)
         XLAL_ERROR(XLAL_EFUNC);
	
	if (approximant == IMRPhenomD)
	{
		double maxamp=0;
		
		REAL8TimeSeries *hp = *hplus;
		REAL8TimeSeries *hc = *hcross;
		int32_t maxind=hp->data->length - 1;
		const REAL8 cfac=cos(inclination);
		const REAL8 pfac = 0.5 * (1. + cfac*cfac);
		for (int32_t loopi =hp->data->length - 1; loopi > -1; loopi--)
		{
			
			 double ampsqr = (hp->data->data[loopi])*(hp->data->data[loopi]) +
					(hc->data->data[loopi])*(hc->data->data[loopi]);
			 if (ampsqr > maxamp)
			 {
					 maxind=loopi;
					 maxamp=ampsqr;
			 }
			 hp->data->data[loopi] *= pfac;
			 hc->data->data[loopi] *= cfac;
		}
	}
  
    /* condition the time domain waveform by tapering in the extra time
     * at the beginning and high-pass filtering above original f_min */
    XLALSimInspiralTDConditionStage1(
		*hplus, 
		*hcross, 
		 extra_time_fraction * chirp_time_upper_bound.seconds, 
		 original_f_min
	);
  
    return 0;
}

void generatePhenomCUDA(
	const Approximant   approximant,
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
		strain[index].x = (float)hplus->data->data[new_waveform_index];
		strain[index].y = (float)hcross->data->data[new_waveform_index];
    }
	
	free(hcross->data->data); free(hplus->data->data);

	*ret_strain = strain;
}

#endif