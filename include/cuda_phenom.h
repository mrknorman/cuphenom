#ifndef _PHENOM_H
#define _PHENOM_H

#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

// Fraction of waveform duration to add as extra time for tapering:
#define EXTRA_TIME_FRACTION 0.1f
// More extra time measured in cycles at the starting frequency:
#define EXTRA_CYCLES 3.0f

#include "phenomd.h"
#include "phenom_functions.h"
#include "phenom_d_functions.h"

void checkTemporalErrors(
    const temporal_properties_s *temporal_properties,
    const int32_t                num_waveforms
    ) {
    
    for (int32_t index = 0; index < num_waveforms; index++) 
    {
        const frequencyUnit_t reference_frequency  = temporal_properties[index].reference_frequency; 
        const frequencyUnit_t frequency_interval   = temporal_properties[index].frequency_interval;    
        const frequencyUnit_t starting_frequency   = temporal_properties[index].starting_frequency;       
        const frequencyUnit_t old_ending_frequency = temporal_properties[index].ending_frequency;   

        // Check inputs for sanity:
        if (reference_frequency.hertz < 0.0f)
        {
            fprintf(
                stderr,
                "%s:\n"
                "Waveform %i: \n"
                "Warning! Original reference frequency (%f) Hz must be positive (or"
                " 0 for \"ignore\" \n",
                __func__,
                index,
                reference_frequency.hertz
            );   
        }
        if (frequency_interval.hertz < 0.0f)
        {
            fprintf(
                stderr,
                "%s:\n"
                "Waveform %i: \n"
                "Warning! Frequency interval (%f) Hz must be positive \n",
                __func__,
                index,
                frequency_interval.hertz
            );   
        }
        if (starting_frequency.hertz <= 0.0f)
        {
            fprintf(
                stderr,
                "%s:\n"
                "Waveform %i: \n"
                "Warning! starting_frequency (%f) Hz must be non-negative.\n",
                __func__,
                index,
                starting_frequency.hertz
            );   
        }
        if (old_ending_frequency.hertz < 0.0f)
        {
            fprintf(
                stderr,
                "%s:\n"
                "Waveform %i: \n"
                "Warning! old_ending_frequency (%f) Hz must be positive. \n",
                __func__,
                index,
                old_ending_frequency.hertz
            );   
        }
    }   
}

m_complex_waveform_axes_s cuPhenomDGenerateFD(
    const system_properties_s   *system_properties,
          temporal_properties_s *temporal_properties,
    const int32_t                num_waveforms
    //NRTidal_version_type      NRTidal_version /**< Version of NRTides; can be one of NRTidal versions or NoNRT_V for the BBH baseline */
    ) {
    
    // Check through temporal properties for last min errors:
    checkTemporalErrors(
        temporal_properties,
        num_waveforms
    );
    
    m_complex_waveform_axes_s waveform_axes_fd = 
        initPhenomDWaveformAxes(
            temporal_properties,
            system_properties,
            num_waveforms
        );
                
    sumPhenomDFrequencies(
        waveform_axes_fd
    );
    
    /*
    cudaFree(waveform_axes_fd.strain.values);
    cudaFree(waveform_axes_fd.strain.num_samples_in_waveform);
    cudaFree(waveform_axes_fd.aproximant_variables_of);
    cudaFree(waveform_axes_fd.system_properties_of);
    cudaFree(waveform_axes_fd.temporal_properties_of);
    cudaFree(waveform_axes_fd.merger_time_for_waveform);
        
    cudaFree(waveform_axes_fd.frequency.values);
    cudaFree(waveform_axes_fd.frequency.interval_of_waveform);
    cudaFree(waveform_axes_fd.frequency.num_samples_in_waveform);
    
    cudaFree(waveform_axes_fd.time.values);
    cudaFree(waveform_axes_fd.time.interval_of_waveform);
    cudaFree(waveform_axes_fd.time.num_samples_in_waveform);
    */
    
    for (int32_t index = 0; index < num_waveforms; index++) 
    {
        const frequencyUnit_t reference_frequency  = temporal_properties[index].reference_frequency; 
        const frequencyUnit_t frequency_interval   = temporal_properties[index].frequency_interval;    
        const frequencyUnit_t starting_frequency   = temporal_properties[index].starting_frequency;       
        const frequencyUnit_t old_ending_frequency = temporal_properties[index].ending_frequency;   
    
        // Update temporal properties: (OLD REPLACE WITH GPU SIDE CALCULATION)
        temporal_properties[index].reference_frequency = 
            initFrequencyHertz(        
                (reference_frequency.hertz == 0.0f) ? 
                starting_frequency.hertz : reference_frequency.hertz
            );     

        temporal_properties[index].ending_frequency = 
            initFrequencyHertz(
                f_CUT/system_properties[index].total_mass.seconds
            );
    }
            
    // Coalesce at merger_time = 0:
    // Shift by overall length in time:  
    
    /*
    waveform_axes_fd.merger_time_for_waveform[0] = \
        initTimeSeconds(0.0f);
    waveform_axes_fd.merger_time_for_waveform[0] = 
        addTimes(
            2, 
            waveform_axes_fd.merger_time_for_waveform[0], 
            initTimeSeconds(-1.0f / frequency_interval.hertz)
        );
    */
            
    return waveform_axes_fd;
}

m_complex_waveform_axes_s cuInspiralFD(
    const system_properties_s   *system_properties,
    const temporal_properties_s *original_temporal_properties,
    const int32_t                num_waveforms,
    const approximant_e          approximant               
    ) {  
    
    // Create new variable temporal properties array:
    const size_t temporal_array_size = 
        sizeof(temporal_properties_s) * (size_t)num_waveforms;
        
    temporal_properties_s *temporal_properties = 
        (temporal_properties_s*) malloc(temporal_array_size);
    
    temporal_properties = 
        memcpy(
            (void*)temporal_properties, 
            (void*)original_temporal_properties, 
            temporal_array_size
        );
        
        
    timeUnit_t      time_intervals          [num_waveforms];
    frequencyUnit_t new_starting_frequencies[num_waveforms];
    frequencyUnit_t old_starting_frequencies[num_waveforms];
    timeUnit_t      time_shifts             [num_waveforms];
    
    
    for (int32_t index = 0; index < num_waveforms; index++) 
    {
        // Unpack variables:
        frequencyUnit_t old_starting_frequency = 
            temporal_properties[index].starting_frequency;
        frequencyUnit_t frequency_interval = 
            temporal_properties[index].frequency_interval; 
        frequencyUnit_t ending_frequency = 
            temporal_properties[index].ending_frequency;
        frequencyUnit_t reference_frequency = 
            temporal_properties[index].reference_frequency;

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
                (int32_t)roundf(ending_frequency.hertz / frequency_interval.hertz);

            // If not a power of two:
            if ((ending_frequency_num_samples & (ending_frequency_num_samples - 1))) 
            { 
                frexpf((float)ending_frequency_num_samples, &num_samples_exp);
                nyquist_frequency = 
                    initFrequencyHertz(
                        ldexpf(1.0f, num_samples_exp)*frequency_interval.hertz
                    );

                fprintf(
                    stderr,
                    "%s: \n,"
                    "_max/frequency_interval = %f/%f = %f is not a power of two: "
                    "changing ending_frequency to %f",
                    __func__,
                    ending_frequency.hertz, 
                    frequency_interval.hertz, 
                    ending_frequency.hertz/frequency_interval.hertz, 
                    nyquist_frequency.hertz
                );
            }
        }
        time_intervals[index] = 
            initTimeSeconds(0.5f / nyquist_frequency.hertz);

        // Generate a FD waveform and condition it by applying tapers at frequencies 
        // between a frequency below the requested old_starting_frequency and 
        // new_starting_frequency; also wind the waveform in phase in case it would 
        // wrap around at the merger time:

        // If the requested low frequency is below the lowest Kerr ISCO
        // frequency then change it to that frequency:
        const frequencyUnit_t kerr_isco_frequency = 
            calculateKerrISCOFrequency(system_properties[index]);

        if (old_starting_frequency.hertz > kerr_isco_frequency.hertz)
        {
             old_starting_frequency.hertz = kerr_isco_frequency.hertz;
        }

        // IMR model: estimate plunge and merger time 
        // sometimes these waveforms have phases that
        // cause them to wrap-around an amount equal to
        // the merger-ringodwn time, so we will undo
        // that here:

        // Upper bound on the chirp time starting at old_starting_frequency:
        timeUnit_t chirp_time_upper_bound = 
            InspiralChirpTimeBound(
                old_starting_frequency, 
                system_properties[index]
            );

        // New lower frequency to start the waveform: add some extra early
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
                system_properties[index]
            );

        // Revise (over-)estimate of chirp from new start frequency:
        chirp_time_upper_bound = 
            InspiralChirpTimeBound(
                new_starting_frequency, 
                system_properties[index]
            );

        // We need a long enough segment to hold a whole chirp with some padding 
        // length of the chirp in samples:
        const timeUnit_t extra_time = 
            initTimeSeconds(EXTRA_CYCLES / old_starting_frequency.hertz);

        const timeUnit_t total_time_upper_bound =
            addTimes(
                4,
                chirp_time_upper_bound,
                temporal_properties[index].merge_time_upper_bound,
                temporal_properties[index].ringdown_time_upper_bound,
                scaleTime(
                    extra_time,
                    2.0f
                )
            );

        const int32_t total_num_samples = 
            calcNextPow2(
                total_time_upper_bound.seconds / time_intervals[index].seconds
            );

        const timeUnit_t total_time_rounded = 
            initTimeSeconds(
                (float)total_num_samples * time_intervals[index].seconds
            );

        // Calculate frequency resolution:
        if (frequency_interval.hertz == 0.0f)
        {
            frequency_interval.hertz = 1.0f / total_time_rounded.seconds;
        }
        else if (frequency_interval.hertz > (1.0f / total_time_rounded.seconds))
        {
            fprintf(
                stderr,
                "%s:\n"
                "Specified frequency interval of %f Hz is too large for a chirp of "
                "duration %f s",
                __func__,
                frequency_interval.hertz,
                total_time_rounded.seconds
            );
        }

        // Update temporal properties:
        temporal_properties[index].frequency_interval  = frequency_interval;
        temporal_properties[index].starting_frequency  = new_starting_frequency;
        
        new_starting_frequencies[index] = new_starting_frequency;
        old_starting_frequencies[index] = old_starting_frequency;
        
        // We want to make sure that this waveform will give something
        // sensible if it is later transformed into the time domain:
        // to avoid the end of the waveform wrapping around to the beginning,
        // we shift waveform backwards in time and compensate for this
        // shift by adjusting the epoch:

        time_shifts[index] = 
            addTimes(
                2,
                temporal_properties[index].merge_time_upper_bound,
                temporal_properties[index].ringdown_time_upper_bound
            );
        
    }
        
    // Generate the waveform in the frequency domain starting at 
    // new_starting_frequency:
    m_complex_waveform_axes_s waveform_axes_fd;
    switch (approximant)
    {
        case D:
        // Call the waveform driver routine:
        waveform_axes_fd = 
            cuPhenomDGenerateFD(
                system_properties,
                temporal_properties,
                num_waveforms
                //4
            );
        
        break;

        case XPHM:
            
            /*
            for (int32_t index = 0; index < num_waveforms; index++)
            {
                if (reference_frequency.hertz == 0.0f)
                {
                    // Default reference frequency is minimum frequency:
                    reference_frequency = new_starting_frequency;
                }
            }
            */

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
                phiRef, f_start, ending_frequency, frequency_interval, f_ref, 
                LALparams
                );
            */

        break;

        default:
            fprintf(
                stderr, 
                "%s:\n"
                "Warning! Approximant not implemented. \n",
                __func__
            );
        break;
    }
    
    printf("Here \n");
    
    for (int32_t index = 0; index < num_waveforms; index++)
    {
        // Assign waveform properties:
        waveform_axes_fd.frequency.interval_of_waveform[index] = temporal_properties[index].frequency_interval;
        waveform_axes_fd.time.interval_of_waveform[index]      = time_intervals[index];
    }
    
    printf("Here After \n");

    
    if (system_properties[0].ascending_node_longitude > 0.0f) 
    {   
        applyPolarization(
            waveform_axes_fd
        );
    }

    taperWaveform(
        waveform_axes_fd,
        new_starting_frequencies,
        old_starting_frequencies
    );

    performTimeShifts(
        waveform_axes_fd,
        time_shifts
    );
    
    free(temporal_properties);
    
    return waveform_axes_fd;
}

m_waveform_axes_s cuInspiralTDFromFD(
          system_properties_s   *system_properties,
          temporal_properties_s *temporal_properties,
    const int32_t                num_waveforms,
    const approximant_e          approximant 
    ) {
         
    // Generate the conditioned waveform in the frequency domain note: redshift 
    // factor has already been applied above set frequency_interval = 0 to get a 
    // small enough resolution:
    temporal_properties[0].frequency_interval = initFrequencyHertz(0.0f);
    m_complex_waveform_axes_s waveform_axes_fd =
        cuInspiralFD(
            system_properties,
            temporal_properties,
            num_waveforms,
            approximant
        );
    
    // We want to make sure that this waveform will give something
    // sensible if it is later transformed into the time domain:
    // to avoid the end of the waveform wrapping around to the beginning,
    // we shift waveform backwards in time and compensate for this
    // shift by adjusting the epoch:
    
    timeUnit_t time_shifts[] = {
        temporal_properties[0].extra_time
    };
    
    performTimeShifts(
        waveform_axes_fd,
        time_shifts
    );
        
    m_waveform_axes_s waveform_axes_td = 
        convertWaveformFDToTD(
            waveform_axes_fd
        );         

    // Compute how long a chirp we should have revised estimate of chirp length 
    // from new start frequency:
    temporal_properties[0].chirp_time_upper_bound =  
        InspiralChirpTimeBound(
            temporal_properties[0].starting_frequency, 
            system_properties[0]
        );   
    
    temporal_properties[0].starting_frequency = 
        InspiralChirpStartFrequencyBound(
            scaleTime(
                temporal_properties[0].chirp_time_upper_bound, 
                (1.0f + EXTRA_TIME_FRACTION)
            ),
            system_properties[0]
        );
    
    // Integer number of samples:
    const timeUnit_t time_shift = 
        initTimeSeconds(
              (roundf(
                    temporal_properties[0].extra_time.seconds
                  / temporal_properties[0].time_interval.seconds) 
            * temporal_properties[0].time_interval.seconds)
        );
    
    // Amount to snip off at the end is num_extra_samples:
    const int32_t num_extra_samples = 
        (int32_t) roundf(
            time_shift.seconds / temporal_properties[0].time_interval.seconds
        );
    
    // Snip off extra time at the end of waveform:
    waveform_axes_td.strain.num_samples_in_waveform[0] -= (float)num_extra_samples;
           
    return waveform_axes_td;
}

m_waveform_axes_s generateInspiral(
    system_properties_s   *system_properties,
    temporal_properties_s *temporal_properties, 
    const int32_t          num_waveforms,
    const approximant_e    approximant               
    ) {
    
    // Hard coded constants:
    
    // Starting frequency is overwritten below, so keep original value:
    //const frequencyUnit_t original_starting_frequency = starting_frequency; 
    
    // SEOBNR flag for spin aligned model version. 1 for SEOBNRv1, 2 for SEOBNRv2
    float polarization = system_properties[0].ascending_node_longitude;
    
    // Save original inclination:
    const angularUnit_t original_inclination = system_properties[0].inclination;
    
    switch (approximant)
    {
        case D:
            // Generate TD waveforms with zero inclincation so that amplitude 
            // can be calculated from hplus and hcross, apply 
            // inclination-dependent factors in function below:      
            system_properties[0].inclination = initAngleRadians(0.0);
        break;

        case XPHM:
            polarization = 0.0f;

        break;

        default:
            printf(
                "%s:\n"
                "Aproximant (%i) not supported! Need = %i \n",
                __func__, (int) approximant, (int)D
                );
        break;
    }
    
    // Generate the waveform in the time domain starting at starting frequency:
    m_waveform_axes_s waveform_axes_td = 
        cuInspiralTDFromFD(
            system_properties,
            temporal_properties,
            num_waveforms,
            approximant
        );
    
    // Set inclination to original:
    system_properties[0].inclination = original_inclination;
    
    switch (approximant)
    {
        case D:                        
            // Apply inclination-dependent factors:
            ;waveform_axes_td = inclinationAdjust(
                waveform_axes_td
            );
            break;

        break;

        case XPHM:
            polarization = 0.0;
        break;

        default:
            printf(
                "%s:\n"
                "Aproximant (%i) not supported! \n",
                __func__, (int) approximant
            );
        break;
    }
    
    // R.C.: here's the reference explaining why we perform this
    // rotation https://dcc.ligo.org/LIGO-G1900275:
    waveform_axes_td = 
        polarisationRotation(
            polarization,
            waveform_axes_td
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
        
    return waveform_axes_td;
}

void generatePhenomCUDA(
    const approximant_e     approximant,
    const massUnit_t        mass_1, 
    const massUnit_t        mass_2, 
    const frequencyUnit_t   sample_rate, 
    const timeUnit_t        duration, 
    const angularUnit_t     inclination, 
    const lengthUnit_t      distance, 
          float2_t       **ret_strain
    ) {
    
    const int32_t num_samples = 
        (int32_t)floor(sample_rate.hertz*duration.seconds);
    
    // Setup companion structures:
    
    spin_t spin_1 = 
    {
        .x = 0.0f,
        .y = 0.0f,
        .z = 0.0f
    };
    spin_t spin_2 = 
    {
        .x = 0.0f,
        .y = 0.0f,
        .z = 0.0f
    };
    companion_s companion_a = 
    {
        .mass              = mass_1,
        .spin              = spin_1,
        .quadrapole_moment = 0.0f,
        .lambda            = 0.0f
    };
    companion_s companion_b = 
    {
        .mass              = mass_2,
        .spin              = spin_2,
        .quadrapole_moment = 0.0f,
        .lambda            = 0.0f
    };
    
    angularUnit_t   reference_orbital_phase  = initAngleRadians(0.0f);
    float           ascending_node_longitude = 100.0f;
    float           eccentricity             = 0.0f;
    float           mean_periastron_anomaly  = 0.0f;
    timeUnit_t      time_interval        = 
        initTimeSeconds(1.0f/sample_rate.hertz);
    frequencyUnit_t starting_frequency  = 
        calcMinimumFrequency(
            mass_1, 
            mass_2, 
            duration
        );
    float           redshift            = 0.0f;
    frequencyUnit_t reference_frequency = initFrequencyHertz(0.0f);
    
    // Init property structures:
    const int32_t num_waveforms = 1;
    
    system_properties_s   system_properties[num_waveforms];
    temporal_properties_s temporal_properties[num_waveforms];
    
    for (int32_t index = 0; index < num_waveforms; index++)
    {
        system_properties[index] =
            initBinarySystem(
                companion_a,
                companion_b,
                distance,
                redshift,
                inclination,
                reference_orbital_phase,
                ascending_node_longitude,
                eccentricity, 
                mean_periastron_anomaly
            );

        temporal_properties[index] =
            initTemporalProperties(
                time_interval,       // <-- Sampling interval (timeUnit_t).
                starting_frequency,  // <-- Starting GW frequency (frequencyUnit_t).
                reference_frequency, // <-- Reference GW frequency (frequencyUnit_t).
                system_properties[index],
                approximant
            );
    }
    
    m_waveform_axes_s waveform_axes_td = 
        generateInspiral(
            system_properties,
            temporal_properties,
            num_waveforms,
            approximant
        );
    
    float2_t *strain = NULL;
    cudaToHost(
        (void**)&waveform_axes_td.strain.values[
        (int32_t)waveform_axes_td.strain.num_samples_in_waveform[0] - num_samples - 1], 
        sizeof(float2_t),
        num_samples,
        (void**) &strain
    );
    cudaFree(waveform_axes_td.strain.values);
        
    if (waveform_axes_td.strain.num_samples_in_waveform[0] < num_samples) 
    {    
        fprintf(
            stderr, 
            "Warning! Cuphenom not generating waveforms of desired num_samples."
            "\n"
        );
    }

    *ret_strain = strain;
}

#endif