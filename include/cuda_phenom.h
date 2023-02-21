#ifndef _PHENOM_H
#define _PHENOM_H

#include <omp.h>

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#define LAL_PI 3.141592653589793238462643383279502884f
#define LAL_MTSUN_SI 4.925491025543575903411922162094833998e-6f
#define LAL_MSUN_SI 1.988409902147041637325262574352366540e30f
#define LAL_MRSUN_SI 1.476625061404649406193430731479084713e3f

// Fraction of waveform duration to add as extra time for tapering:
#define EXTRA_TIME_FRACTION 0.1f
// More extra time measured in cycles at the starting frequency:
#define EXTRA_CYCLES 3.0f

#include "phenomd.h"
#include "phenom_functions.h"

const LALUnit lalStrainUnit = { 0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALUnit lalSecondUnit = {0,{0,0,1,0,0,0,0},{0,0,0,0,0,0,0}};

// If x and X are approximately equal to relative accuracy epsilon
// then set x = X.
// If X = 0 then use an absolute comparison.

float Nudge(
    const float x, 
    const float X, 
    const float epsilon
    ) {
    
    float new_x = x;
    if (X != 0.0)
    {
        if (!gsl_fcmp(x, X, epsilon))
        {
            printf("Nudging value %.15g to %.15g.\n", x, X);
            new_x = X;
        }
    }
    else
    {
        if (fabs(x - X) < epsilon) 
        {
            new_x = X;
        }
    }
    
    return new_x;
}

complex_waveform_axes_s _cuPhenomDGenerateFD(
     const system_properties_s system_properties,
     const int32_t             num_strain_axis_samples,
     const frequencyUnit_t     starting_frequency,
     const frequencyUnit_t     ending_frequency,
     const frequencyUnit_t     frequency_interval,      
     const frequencyUnit_t     reference_frequency,  // reference frequency 
     NRTidal_version_type      NRTidal_version               // NRTidal version; either NRTidal_V or NRTidalv2_V or NoNRT_V in case of BBH baseline
) {
    
    /*
     * Internal cuPhenom function to generate IMRPhenomD waveform in the
     * frequency domain.
     * @param 
     * @see 
     * @return void
     */
        
    const float non_gdr_chi_6 = 0.0f;
    const int32_t tidal_pn_order = -1; 
    
    // Unpack system_properties:
    const massUnit_t total_mass = system_properties.total_mass;
    
    const float chi1 = system_properties.companion[0].spin.z;
    const float chi2 = system_properties.companion[1].spin.z;
    
    const float distance    = system_properties.distance.meters;
    const float inclination = system_properties.inclination.radians;
    // Phase at reference frequency
    const angularUnit_t reference_phase = 
        system_properties.reference_orbital_phase;    

            
    // Check frequency bounds:
    if (starting_frequency.hertz <= 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n Error! Starting frequency (%f) Hz must be positive. \n",
            __func__, starting_frequency.hertz 
        );
    }
    if (ending_frequency.hertz < 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n Error! Ending frequency (%f) Hz must not be negative. \n",
            __func__, ending_frequency.hertz 
        );
    }
    
    // Init symmetric mass ratio:
    float symmetric_mass_ratio = system_properties.symmetric_mass_ratio;
    
    if (symmetric_mass_ratio > 0.25f)
    {
        symmetric_mass_ratio = Nudge(symmetric_mass_ratio, 0.25f, 1.0e-6f);
    }
    if (symmetric_mass_ratio > 0.25f || symmetric_mass_ratio < 0.0f)
    {
        fprintf(
            stderr, 
            "%s:\n"
            "Unphysical symmetric_mass_ratio. Must be between 0. and 0.25.\n", 
            __func__
        );
    }

    // Compute the amplitude pre-factor:
    const float amp0 = 
          2.0f 
        * sqrtf(5.0f / (64.0f*(float)M_PI)) 
        * total_mass.meters 
        * total_mass.seconds/distance;
    
    // Calculate index shift between freqs and the frequency series:
    int32_t offset = 0; 
    if (frequency_interval.hertz > 0.0f) 
    { 
        offset = (int32_t) (starting_frequency.hertz / frequency_interval.hertz);
    }
    
    // Initilise waveform axes:
    complex_waveform_axes_s waveform_axes =
        generatePhenomD(
            starting_frequency,
            ending_frequency,
            frequency_interval,
            num_strain_axis_samples
        );
    
    // Calculate phenomenological parameters:
    const float final_spin = calcFinalSpin(system_properties); 

    if (final_spin < MIN_FINAL_SPIN)
    {
        fprintf(
            stderr, 
            "%s: \n"
            "Final spin (%f) and ISCO frequency of this system_properties"
            " are small, the model might misbehave here.",
            __func__, final_spin
        );
    }
    
    // Init amplitude coefficients:
    amplitude_coefficients_s amplitude_coefficients = 
        initAmplitudeCoefficients(
            symmetric_mass_ratio, 
            chi1, 
            chi2, 
            final_spin
        );
    
    phase_coefficients_s phase_coefficients = 
        initPhaseCoefficients(
            symmetric_mass_ratio, 
            chi1, 
            chi2, 
            final_spin
        );
    
    pn_phasing_series_s phasing_series =
        initTaylorF2AlignedPhasingSeries(system_properties, tidal_pn_order);
    
    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 
    // implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but was 
    // not available when PhenomD was tuned.
    const float testGRcor = 1.0f + non_gdr_chi_6;
    
    phasing_series.v[6] -= 
        subtract3PNSS(system_properties)*phasing_series.v[0]*testGRcor;
    
    phase_inspiral_prefactors_s phase_prefactors = 
        initPhaseInspiralPrefactors(phase_coefficients, phasing_series);
    
    const float Rholm = 1.0f;
    const float Taulm = 1.0f;
    
    // Compute coefficients to make phase C^1 continuous 
    // (phase and first derivative):
    phase_coefficients = 
        computePhenomDPhaseConnectionCoefficients(
            phase_coefficients, 
            phasing_series, 
            phase_prefactors, 
            Rholm, 
            Taulm
        );

    //time shift so that peak amplitude is approximately at t=0
    //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
    const float phase_shift = 
        calculateMergerRingdownPhaseAnsatzDerivitive(
            amplitude_coefficients.inspiral_merger_peak_frequency, 
            phase_coefficients, 
            Rholm, 
            Taulm
        );

    const amplitude_inspiral_prefactors_s amplitude_prefactors = 
        initAmplitudeInspiralPrefactors(amplitude_coefficients);

    // Incorporating reference frequency:
    const float reference_mass_frequency = 
        total_mass.seconds * reference_frequency.hertz;
    useful_powers_s powers_of_reference_mass =
        initUsefulPowers(reference_mass_frequency);
    
    const float phase = 
        calculatePhase(
            powers_of_reference_mass, 
            phase_coefficients, 
            phase_prefactors,
            Rholm, 
            Taulm
        );

    // Factor of 2 b/c reference_phase is orbital phase:
    const float phi_precalc = 2.0f*reference_phase.radians + phase;
    
    sumPhenomDFrequencies(
        waveform_axes,
        inclination,
        total_mass.seconds,
        amplitude_coefficients,
        amplitude_prefactors,
        phase_coefficients, 
        phase_prefactors, 
        offset,
        phase_shift,
        amp0,
        reference_mass_frequency,
        phi_precalc
    );
    
    cudaFree(waveform_axes.frequency.values);        
    return waveform_axes;
}

complex_waveform_axes_s cuPhenomDGenerateFD(
    const system_properties_s system_properties,
    const frequencyUnit_t     original_reference_frequency, /**< reference frequency (Hz) */
    const frequencyUnit_t     frequency_interval,                /**< Sampling frequency (Hz) */
    const frequencyUnit_t     starting_frequency,                 /**< Starting GW frequency (Hz) */
    const frequencyUnit_t     old_ending_frequency,                 /**< End frequency; 0 defaults to Mf = \ref f_CUT */
    NRTidal_version_type      NRTidal_version /**< Version of NRTides; can be one of NRTidal versions or NoNRT_V for the BBH baseline */
    ) {
    
    // Check inputs for sanity:
    if (original_reference_frequency.hertz < 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! Original reference frequency (%f) Hz must be positive (or"
            " 0 for \"ignore\" \n",
            __func__,
            original_reference_frequency.hertz
        );   
    }
    if (frequency_interval.hertz < 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! Frequency interval (%f) Hz must be positive \n",
            __func__,
            frequency_interval.hertz
        );   
    }
    if (starting_frequency.hertz <= 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! starting_frequency (%f) Hz must be non-negative.\n",
            __func__,
            starting_frequency.hertz
        );   
    }
    if (old_ending_frequency.hertz < 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! old_ending_frequency (%f) Hz must be positive. \n",
            __func__,
            old_ending_frequency.hertz
        );   
    }
           
    // If no reference frequency given, set it to the starting GW frequency:
    const frequencyUnit_t new_reference_frequency = 
        initFrequencyHertz(        
            (original_reference_frequency.hertz == 0.0f) ? 
            starting_frequency.hertz : original_reference_frequency.hertz
        );

    const float frequency_cutoff = f_CUT/system_properties.total_mass.seconds;
        
    // Somewhat arbitrary end point for the waveform.
    // Chosen so that the end of the waveform is well after the ringdown.
    if (frequency_cutoff <= starting_frequency.hertz) 
    {
        fprintf(
           stderr,
           "%s:\n"
           "Warning! (frequency_cutoff: %fHz) <= (starting_frequency = %fHz) \n. \n",
            __func__,
            frequency_cutoff,
            starting_frequency.hertz
        );
    }

    // Default old_ending_frequency to cuttoff_frequency:
    float new_ending_frequency_hertz = old_ending_frequency.hertz;
    new_ending_frequency_hertz = old_ending_frequency.hertz ? old_ending_frequency.hertz : frequency_cutoff;
    new_ending_frequency_hertz = 
        (new_ending_frequency_hertz > frequency_cutoff) ? frequency_cutoff : new_ending_frequency_hertz;
    
    if (new_ending_frequency_hertz <= starting_frequency.hertz)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! new_ending_frequency_hertz (%f) <= starting_frequency (%f)\n",
            __func__,
            new_ending_frequency_hertz,
            starting_frequency.hertz
        );
    }
    
    const frequencyUnit_t new_ending_frequency = 
        initFrequencyHertz(new_ending_frequency_hertz);
    
    int32_t num_strain_axis_samples = 0; 
    timeUnit_t gps_time = initTimeSeconds(0.0f);

    // Set up output array with size closest power of 2:
    num_strain_axis_samples = calcNextPow2(old_ending_frequency.hertz / frequency_interval.hertz) + 1;

    // Coalesce at gps_time = 0:
    // Shift by overall length in time:  
    gps_time = addTimes(2, gps_time, initTimeSeconds(-1.0f / frequency_interval.hertz));
    
    if (new_ending_frequency.hertz < old_ending_frequency.hertz) 
    {
        // The user has requested a higher old_ending_frequency than Mf=frequency_cutoff.resize 
        // the frequency series to fill with zeros beyond the cutoff frequency:
        num_strain_axis_samples = calcNextPow2(old_ending_frequency.hertz / frequency_interval.hertz) + 1; 
        // We actually want to have the length be a power of 2 + 1.
    }
        
    complex_waveform_axes_s waveform_axes = 
        _cuPhenomDGenerateFD(
            system_properties,
            num_strain_axis_samples,
            starting_frequency,
            new_ending_frequency,
            frequency_interval, 
            new_reference_frequency,
            NRTidal_version
        );
    
    return waveform_axes;
}

void _convertWaveformAxesToLAL_temp(
    complex_waveform_axes_s waveform_axes,
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde     /**< FD cross polarization */
    ) {
    
    // Produce both polarizations:
    const LALUnit lalStrainUnit = { 0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
    
    // Freqs contains uniform frequency grid with spacing frequency_interval; we 
    // start at frequency 0:
    if (waveform_axes.frequency.interval.hertz > 0.0f) 
    { 
        *hptilde = 
            XLALCreateCOMPLEX16FrequencySeries(
                "htilde: FD waveform", 
                (LIGOTimeGPS){ 0, 0 }, 
                0.0, 
                waveform_axes.frequency.interval.hertz, 
                &lalStrainUnit, 
                waveform_axes.strain.num_samples
            );

        *hctilde = 
            XLALCreateCOMPLEX16FrequencySeries(
                "ctilde: FD waveform", 
                (LIGOTimeGPS){ 0, 0 }, 
                0.0, 
                waveform_axes.frequency.interval.hertz, 
                &lalStrainUnit, 
                waveform_axes.strain.num_samples
            );
    }
    else // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
    {
        *hptilde = (COMPLEX16FrequencySeries*)
            XLALCreateCOMPLEX16FrequencySeries(
                "htilde: FD waveform", 
                (LIGOTimeGPS){ 0, 0 },
                waveform_axes.frequency.values[0].hertz, 
                waveform_axes.frequency.interval.hertz, 
                &lalStrainUnit, 
                waveform_axes.strain.num_samples
            );
        *hctilde = (COMPLEX16FrequencySeries*)
            XLALCreateCOMPLEX16FrequencySeries(
                "ctilde: FD waveform", 
                (LIGOTimeGPS){ 0, 0 }, 
                waveform_axes.frequency.values[0].hertz, 
                waveform_axes.frequency.interval.hertz, 
                &lalStrainUnit, 
                waveform_axes.strain.num_samples
            );
    }

    complex_strain_element_t *strain_fd = NULL;
    cudaToHost(
        (void*)waveform_axes.strain.values, 
        sizeof(complex_strain_element_t),
        waveform_axes.strain.num_samples,
        (void**)&strain_fd
    );
    
    for(int32_t j = 0; j < (int32_t)(*hptilde)->data->length; j++) 
    {    
        (*hctilde)->data->data[j] = (complex double)strain_fd[j].cross;
        (*hptilde)->data->data[j] = (complex double)strain_fd[j].plus;
    }
}

complex_waveform_axes_s cuInspiralFD(
    const system_properties_s system_properties,
    frequencyUnit_t frequency_interval,     // < sampling interval (Hz)
    frequencyUnit_t old_starting_frequency, // < starting GW frequency (Hz)
    frequencyUnit_t ending_frequency,       // < ending GW frequency (Hz)
    frequencyUnit_t reference_frequency,    // < Reference frequency (Hz)
    Approximant approximant                 // < post-Newtonian approximant to use for waveform production
    ) {  
    
    // Apply condition that ending_frequency rounds to the next power-of-two 
    // multiple of frequency_interval.
    // Round ending_frequency / frequency_interval to next power of two.
    // Set ending_frequency to the new Nyquist frequency.
    // The length of the chirp signal is then 2 * nyquist_frequency / 
    // frequency_interval.
    // The time spacing is 1 / (2 * nyquist_frequency):
    frequencyUnit_t nyquist_frequency = ending_frequency;
    int32_t ending_frequency_num_samples = 0, num_samples_exp = 0;
    if (frequency_interval.hertz != 0) 
    {
        ending_frequency_num_samples = 
            (int32_t)round(ending_frequency.hertz / frequency_interval.hertz);
                
        // If not a power of two:
        if ((ending_frequency_num_samples & (ending_frequency_num_samples - 1))) 
        { 
            frexpf(ending_frequency_num_samples, &num_samples_exp);
            nyquist_frequency = 
                initFrequencyHertz(
                    ldexpf(1.0f, num_samples_exp)*frequency_interval.hertz
                );
            
            fprintf(
                stderr,
                "%s: \n,"
                "_max/frequency_interval = %f/%f = %f is not a power of two: "
                "changing ending_frequency to %f",
                ending_frequency.hertz, 
                frequency_interval.hertz, 
                ending_frequency.hertz/frequency_interval.hertz, 
                nyquist_frequency.hertz
            );
        }
    }
    timeUnit_t sampling_interval = initTimeSeconds(0.5f / nyquist_frequency.hertz);

    // Generate a FD waveform and condition it by applying tapers at frequencies 
    // between a frequency below the requested old_starting_frequency and 
    // new_starting_frequency; also wind the waveform in phase in case it would 
    // wrap around at the merger time:
    
    // If the requested low frequency is below the lowest Kerr ISCO
    // frequency then change it to that frequency:
    const frequencyUnit_t kerr_isco_frequency = 
        calculateKerrISCOFrequency(system_properties);
    
    if (old_starting_frequency.hertz > kerr_isco_frequency.hertz)
    {
         old_starting_frequency.hertz = kerr_isco_frequency.hertz;
    }
    
    // IMR model: estimate plunge and merger time 
    // sometimes these waveforms have phases that
    // cause them to wrap-around an amount equal to
    // the merger-ringodwn time, so we will undo
    // that here:
    const timeUnit_t merge_time_upper_bound = 
        addTimes(
            2,
            InspiralMergeTimeBound(system_properties),
            InspiralRingdownTimeBound(system_properties)
        );
    
    // Upper bound on the chirp time starting at old_starting_frequency:
    timeUnit_t chirp_time_upper_bound = 
        InspiralChirpTimeBound(
            old_starting_frequency, 
            system_properties
        );
    
    //  new lower frequency to start the waveform: add some extra early
    // part over which tapers may be applied, the extra amount being
    // a fixed fraction of the chirp time; add some additional padding
    // equal to a few extra cycles at the low frequency as well for
    // safety and for other routines to use: 
    const frequencyUnit_t new_starting_frequency = 
        InspiralChirpStartFrequencyBound(
            scaleTime(
                chirp_time_upper_bound, 
                (1.0f + EXTRA_TIME_FRACTION)
            ),
            system_properties
        );
    
    // Revise (over-)estimate of chirp from new start frequency:
    chirp_time_upper_bound = 
        InspiralChirpTimeBound(
            new_starting_frequency, 
            system_properties
        );

    // We need a long enough segment to hold a whole chirp with some padding 
    // length of the chirp in samples:
    const timeUnit_t extra_time = 
        initTimeSeconds(EXTRA_CYCLES / old_starting_frequency.hertz);
    
    const timeUnit_t total_time_upper_bound =
        addTimes(
            3,
            chirp_time_upper_bound,
            merge_time_upper_bound,
            scaleTime(
                extra_time,
                2.0f
            )
        );
    
    const float total_num_samples = 
        calcNextPow2(total_time_upper_bound.seconds / sampling_interval.seconds);
    
    // Frequency resolution:
    if (frequency_interval.hertz == 0.0f)
    {
        frequency_interval.hertz = 1.0f / (total_num_samples * sampling_interval.seconds);
    }
    else if (frequency_interval.hertz > 1.0f / (total_num_samples * sampling_interval.seconds))
    {
        fprintf(
            stderr,
            "%s:\n"
            "Specified frequency interval of %f Hz is too large for a chirp of "
            "duration %f s",
            frequency_interval.hertz,
            total_num_samples * sampling_interval.seconds
        );
    }
    
    // Generate the waveform in the frequency domain starting at 
    // new_starting_frequency:
    
    complex_waveform_axes_s waveform_axes;
    switch (approximant)
    {
        case IMRPhenomD:
        // Call the waveform driver routine:
        waveform_axes = 
            cuPhenomDGenerateFD(
                system_properties,
                reference_frequency, 
                frequency_interval, 
                new_starting_frequency, 
                ending_frequency, 
                4
            );
        
        break;

        case IMRPhenomXPHM:
            if (reference_frequency.hertz == 0.0f)
            {
                // Default reference frequency is minimum frequency:
                reference_frequency = new_starting_frequency;
            }

            // Call the main waveform driver. Note that we pass the full spin 
            // vectors with XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame 
            // being effectively called in the initialization of the pPrec 
            // struct
            
            /* Disabled for now
            XLALSimIMRPhenomXPHM(
                hptilde, hctilde,
                m1, m2,
                S1x, S1y, S1z,
                S2x, S2y, S2z,
                distance, inclination,
                phiRef, f_start, ending_frequency, frequency_interval, f_ref, LALparams
                );
            */

        break;

        default:
            fprintf(
                stderr, 
                "FD version of approximant not implemented in lalsimulation\n"
            );
        break;
    }
    
    waveform_axes.frequency.interval = frequency_interval;

    if (system_properties.ascending_node_longitude > 0.0f) 
    {   
        applyPolarization(
            waveform_axes,
            system_properties.ascending_node_longitude
        );
    }
        
    taperWaveform(
        waveform_axes,
        new_starting_frequency.hertz,
        old_starting_frequency.hertz,
        frequency_interval.hertz
    );
    
    // Integer number of time samples:
    const float num_samples_time_shift =
          roundf(merge_time_upper_bound.seconds / sampling_interval.seconds) 
        * sampling_interval.seconds; 
    
    // We want to make sure that this waveform will give something
    // sensible if it is later transformed into the time domain:
    // to avoid the end of the waveform wrapping around to the beginning,
    // we shift waveform backwards in time and compensate for this
    // shift by adjusting the epoch:
    performTimeShift(
          waveform_axes,
          num_samples_time_shift
    );

    //XLALGPSAdd(&(*hptilde)->epoch, tshift);
    //XLALGPSAdd(&(*hctilde)->epoch, tshift);

    return waveform_axes;
}

int cuInspiralTDFromFD(
     REAL8TimeSeries       **hplus,                    /**< +-polarization waveform */
     REAL8TimeSeries       **hcross,                   /**< x-polarization waveform */
     system_properties_s     system_properties,
     temporal_properties_s   temporal_properties,
     LALDict *LALparams,                         /**< LAL dictionary containing accessory parameters */
     Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
) {
    
        
    // Generate the conditioned waveform in the frequency domain note: redshift 
    // factor has already been applied above set deltaF = 0 to get a small
    // enough resolution:
    complex_waveform_axes_s waveform_axes = 
        cuInspiralFD(
            system_properties,
            initFrequencyHertz(0.0f), 
            temporal_properties.starting_frequency, 
            temporal_properties.ending_frequency, 
            temporal_properties.reference_frequency,
            approximant
        );
    
    float num_samples_time_shift = 
        roundf(temporal_properties.extra_time.seconds/temporal_properties.sampling_interval.seconds) 
        * temporal_properties.sampling_interval.seconds; // integer number of samples 
    
    // We want to make sure that this waveform will give something
    // sensible if it is later transformed into the time domain:
    // to avoid the end of the waveform wrapping around to the beginning,
    // we shift waveform backwards in time and compensate for this
    // shift by adjusting the epoch:
    performTimeShift(
          waveform_axes,
          num_samples_time_shift
    );
    
    COMPLEX16FrequencySeries *hptilde = NULL;
    COMPLEX16FrequencySeries *hctilde = NULL;
    
    _convertWaveformAxesToLAL_temp(
        waveform_axes,
        &hptilde,     //< FD plus polarization 
        &hctilde     //< FD cross polarization 
    );
    
    double tshift;

    //XLALGPSAdd(&hptilde->epoch, tshift);
    //XLALGPSAdd(&hctilde->epoch, tshift);
    
    const LALUnit lalStrainUnit = { 0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
    // Rransform the waveform into the time domain:
    size_t chirplen = 2 * (hptilde->data->length - 1);
    *hplus = XLALCreateREAL8TimeSeries("H_PLUS", &hptilde->epoch, 0.0, temporal_properties.sampling_interval.seconds, &lalStrainUnit, chirplen);
    *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &hctilde->epoch, 0.0, temporal_properties.sampling_interval.seconds, &lalStrainUnit, chirplen);
    
    size_t num_waveform_samples = 2 * (hptilde->data->length - 1);

    performIRFFT64(
        hptilde->data->data,
        hctilde->data->data,
        (*hplus)->data->data,
        (*hcross)->data->data,
        temporal_properties.sampling_interval,
        (int32_t)num_waveform_samples
    );
        
    /*
        performIRFFT(
            hptilde->data->data,
            hctilde->data->data,
            (*hplus)->data->data,
            (*hcross)->data->data,
            temporal_properties,
            initFrequencyHertz(hptilde->deltaF),
            (int32_t)num_waveform_samples
        );
    */
    
    
     /*
     XLALGPSAdd(&hptilde->epoch, tshift);
     XLALGPSAdd(&hctilde->epoch, tshift);
  
     // transform the waveform into the time domain 
     size_t chirplen = 2 * (hptilde->data->length - 1);
     *hplus = XLALCreateREAL8TimeSeries("H_PLUS", &hptilde->epoch, 0.0, temporal_properties.sampling_interval.seconds, &lalStrainUnit, chirplen);
     *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &hctilde->epoch, 0.0, temporal_properties.sampling_interval.seconds, &lalStrainUnit, chirplen);
     void* plan = XLALCreateReverseREAL8FFTPlan(chirplen, 0);
     if (!(*hplus) || !(*hcross) || !plan) {
         XLALDestroyCOMPLEX16FrequencySeries(hptilde);
         XLALDestroyCOMPLEX16FrequencySeries(hctilde);
         XLALDestroyREAL8TimeSeries(*hcross);
         XLALDestroyREAL8TimeSeries(*hplus);
         XLALDestroyREAL8FFTPlan(plan);
         //XLAL_ERROR(XLAL_EFUNC);
     }
     XLALREAL8FreqTimeFFT(*hplus, hptilde, plan);
     XLALREAL8FreqTimeFFT(*hcross, hctilde, plan);
    */
    
    /* compute how long a chirp we should have */
    /* revised estimate of chirp length from new start frequency */
    
    temporal_properties.chirp_time_upper_bound =  
        InspiralChirpTimeBound(
            temporal_properties.starting_frequency, 
            system_properties
        );   
    
    temporal_properties.starting_frequency = 
        InspiralChirpStartFrequencyBound(
            scaleTime(
                temporal_properties.chirp_time_upper_bound, 
                (1.0f + EXTRA_TIME_FRACTION)
            ),
            system_properties
        );
    
    // WAT??
    const timeUnit_t new_inspiral_time_upper_bound = 
        InspiralChirpTimeBound(
            temporal_properties.starting_frequency, 
            system_properties
        );

    /* total expected chirp length includes merger */
    chirplen = 
        (size_t) roundf(
            (new_inspiral_time_upper_bound.seconds 
            + temporal_properties.chirp_time_upper_bound.seconds)
            / temporal_properties.sampling_interval.seconds
        );
    
    tshift = (double)
          (roundf(temporal_properties.extra_time.seconds / temporal_properties.sampling_interval.seconds) 
        * temporal_properties.sampling_interval.seconds); // Integer number of samples
    
    /* amount to snip off at the end is tshift */
    size_t end = (*hplus)->data->length - (size_t)round(tshift / (double)temporal_properties.sampling_interval.seconds);

    /* snip off extra time at beginning and at the end */
    XLALResizeREAL8TimeSeries(*hplus,  (int32_t) (end - chirplen), chirplen);
    XLALResizeREAL8TimeSeries(*hcross, (int32_t) (end - chirplen), chirplen);

    /* clean up */
    //XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    //XLALDestroyCOMPLEX16FrequencySeries(hctilde);
     
    return 0;
}

void polarisationRotation(
    const double  polarization,
          double *hplus,
          double *hcross,
    const int32_t num_samples
    ) {
    
    //R.C.: here's the reference explaining why we perform this rotation https://dcc.ligo.org/LIGO-G1900275
    if (polarization != 0.0) 
    {        
        double temp_plus = 0.0;
        double cos_polarization = cos(2.0*polarization);
        double sin_polarization = sin(2.0*polarization);
        
        for (int32_t index = 0; index < num_samples; index++) 
        {
            temp_plus = hplus[index];
            hplus [index] = (cos_polarization*hplus [index])
                          + (sin_polarization*hcross[index]);
                          
            hcross[index] = (cos_polarization*hcross[index]) 
                          - (sin_polarization*temp_plus);
        }
    }
}

int32_t inclinationAdjust(
          system_properties_s  system_properties,
          double              *hplus,
          double              *hcross,
    const int32_t              num_samples
    ) {
    
    double  max_amplitude = 0.0;
    int32_t max_index     = num_samples - 1;
    
    const double cross_factor = cos(system_properties.inclination.radians);
    const double plus_factor  = 0.5 * (1.0 + cross_factor*cross_factor);

    for (int32_t index = num_samples - 1; index > -1; index--)
    {
        double amplitude_squared = 
            hplus[index]*hplus[index] + hcross[index]*hcross[index];
            
        if (amplitude_squared > max_amplitude)
        {
            max_index     = index;
            max_amplitude = amplitude_squared;
        }
        
        hplus [index] *= plus_factor;
        hcross[index] *= cross_factor;
    }
    
    return max_index;
}

int32_t generateInspiral(
    // +-polarization waveform:
    REAL8TimeSeries **hplus, 
    // x-polarization waveform:
    REAL8TimeSeries **hcross,    
    // Structure containing properties of companion a:
    companion_s       companion_a,     
    // Structure containing properties of companion b:
    companion_s       companion_b,    
    // Distance of source (lengthUnit_t):
    lengthUnit_t      distance,      
    // Redshift of source:
    double            redshift,       
    // Inclination of source (angularUnit_t):
    angularUnit_t     inclination,         
    // Reference orbital phase (angularUnit_t):
    angularUnit_t     reference_orbital_phase, 
    // longitude of ascending nodes, degenerate with the polarization angle:
    double            ascending_node_longitude,
    // Eccentrocity at reference epoch:
    double            eccentricity,           
    // Mean anomaly of periastron:
    double            mean_periastron_anomaly,  
    // Sampling interval (timeUnit_t):
    timeUnit_t        sampling_interval,     
    // Starting GW frequency (frequencyUnit_t):
    frequencyUnit_t   starting_frequency,    
    // Reference GW frequency (frequencyUnit_t):
    frequencyUnit_t   reference_frequency,     
    // LAL dictionary containing accessory parameters:
    LALDict          *LALparams,                
    // Post-Newtonian approximant to use for waveform production:
    Approximant       approximant               
    ) {
    
    // Hard coded constants:
    
    // Starting frequency is overwritten below, so keep original value:
    //const frequencyUnit_t original_starting_frequency = starting_frequency; 
    
    // Init property structures:
    system_properties_s system_properties =
        initBinarySystem(
            companion_a,
            companion_b,
            distance,
            (float) redshift,
            inclination,
            reference_orbital_phase,
            (float) ascending_node_longitude,
            (float) eccentricity, 
            (float) mean_periastron_anomaly
        );
    
    temporal_properties_s temporal_properties =
        initTemporalProperties(
            sampling_interval,   // <-- Sampling interval (timeUnit_t).
            starting_frequency,  // <-- Starting GW frequency (frequencyUnit_t).
            reference_frequency, // <-- Reference GW frequency (frequencyUnit_t).
            system_properties,
            approximant
        );
    
    // SEOBNR flag for spin aligned model version. 1 for SEOBNRv1, 2 for SEOBNRv2
    double polarization = system_properties.ascending_node_longitude;
    
    // Save original inclination:
    const angularUnit_t original_inclination = system_properties.inclination;
    
    switch (approximant)
    {
        case IMRPhenomD:
            // Generate TD waveforms with zero inclincation so that amplitude 
            // can be calculated from hplus and hcross, apply 
            // inclination-dependent factors in function below:      
            system_properties.inclination = initAngleRadians(0.0);
        break;

        case IMRPhenomXPHM:
            polarization = 0.0;

        break;

        default:
            printf("Aproximant not supported! \n");
        break;
    }
    
    // Generate the waveform in the time domain starting at starting frequency:
    int32_t return_value = cuInspiralTDFromFD(
        hplus, 
        hcross, 
        system_properties,
        temporal_properties,
        LALparams, 
        approximant
    );
    
    // Set inclination to original:
    system_properties.inclination = original_inclination;
    
    switch (approximant)
    {
        case IMRPhenomD:                        
            // Apply inclination-dependent factors:
            ;const int32_t max_index = 
            inclinationAdjust(
                system_properties,
                (*hplus)->data->data,
                (*hcross)->data->data,
                (int32_t) (*hplus)->data->length
            );

            //XLALGPSSetREAL8(&(hp->epoch), -1.0*temporal_properties.sampling_interval.seconds * max_index);
            //XLALGPSSetREAL8(&(hc->epoch), -1.0*temporal_properties.sampling_interval.seconds * max_index);    
        break;

        case IMRPhenomXPHM:
            polarization = 0.0;
        break;

        default:
            printf("Aproximant not supported! \n");
        break;
    }
    
    // R.C.: here's the reference explaining why we perform this
    // rotation https://dcc.ligo.org/LIGO-G1900275:
    polarisationRotation(
        polarization,
        (*hplus)->data->data,
        (*hcross)->data->data,
        (int32_t) (*hplus)->data->length
    );
    
    // Condition the time domain waveform by tapering in the extra time
    // at the beginning and high-pass filtering above original 
    // starting_frequency:
    
    //Leave out for now
    /*
    XLALSimInspiralTDConditionStage1(
        *hplus, 
        *hcross, 
        EXTRA_TIME_FRACTION*temporal_properties.chirp_time_upper_bound.seconds
        + temporal_properties.extra_time.seconds, 
        original_starting_frequency.hertz
    );
    */
        
    return return_value;
}

void generatePhenomCUDA(
    const Approximant       approximant,
    const massUnit_t        mass_1, 
    const massUnit_t        mass_2, 
    const frequencyUnit_t   sample_rate, 
    const timeUnit_t        duration, 
    const angularUnit_t     inclination, 
    const lengthUnit_t      distance, 
          float2_t     **ret_strain
    ) {
    
    const int32_t num_samples = 
        (int32_t)floor(sample_rate.hertz*duration.seconds);
    
    REAL8TimeSeries *hplus  = NULL;
    REAL8TimeSeries *hcross = NULL;
    
    // Setup companion structures:
    
    spin_t spin_1 = 
    {
        .x = 0.0,
        .y = 0.0,
        .z = 0.0
    };
    spin_t spin_2 = 
    {
        .x = 0.0,
        .y = 0.0,
        .z = 0.0
    };
    companion_s companion_1 = 
    {
        .mass              = mass_1,
        .spin              = spin_1,
        .quadrapole_moment = 0.0,
        .lambda            = 0.0
    };
    companion_s companion_2 = 
    {
        .mass              = mass_2,
        .spin              = spin_2,
        .quadrapole_moment = 0.0,
        .lambda            = 0.0
    };
    
    angularUnit_t   reference_orbital_phase  = initAngleRadians(0.0);
    double          ascending_node_longitude = 100.0;
    double          eccentricity             = 0.0;
    double          mean_periastron_anomaly  = 0.0;
    timeUnit_t      sampling_interval        = 
        initTimeSeconds(1.0f/sample_rate.hertz);
    frequencyUnit_t starting_frequency  = 
        calcMinimumFrequency(
            mass_1, 
            mass_2, 
            duration
        );
    
    double          redshift            = 0.0;
    frequencyUnit_t reference_frequency = initFrequencyHertz(0.0);
    LALDict         *extraParams        = NULL;
    
    generateInspiral(
        &hplus,
        &hcross,
        companion_1,
        companion_2,
        distance,
        redshift,
        inclination,
        reference_orbital_phase,
        ascending_node_longitude,
        eccentricity,
        mean_periastron_anomaly,
        sampling_interval,
        starting_frequency,
        reference_frequency,
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
    
    size_t new_array_size = (size_t)num_samples * sizeof(float2_t);

    float2_t *strain = (float2_t*)malloc(new_array_size);
    int32_t new_waveform_index = 0;
    for (int32_t index = 0; index < num_samples; index++) 
    {    
        new_waveform_index = waveform_num_samples - num_samples - 1 + index;
        strain[index].x = (float) hplus->data->data[new_waveform_index];
        strain[index].y = (float) hcross->data->data[new_waveform_index];    
    }
    
    free(hcross->data->data); free(hplus->data->data);

    *ret_strain = strain;
}

#endif