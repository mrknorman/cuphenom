#ifndef _PyCuPhenom_H
#define _PyCuPhenom_H

#include "cuda_phenom.h"

float *pythonWrapperPhenomD(
    const int    approximant_enum,
    const float  mass_1_msun,
    const float  mass_2_msun,
    const float  sample_rate_hertz,
    const float  duration_seconds,
    const float  inclination_radians,
    const float  distance_mpc,
    const float  reference_orbital_phase_in,
    const float  ascending_node_longitude,
    const float  eccentricity,
    const float  mean_periastron_anomaly,
    const float *spin_1_in,
    const float *spin_2_in
    ) {
    
    approximant_e approximant = D; 
    
    const massUnit_t      mass_1       = initMassSolarMass(mass_1_msun);
    const massUnit_t      mass_2       = initMassSolarMass(mass_2_msun);
    const frequencyUnit_t sample_rate  = initFrequencyHertz(sample_rate_hertz);
    const timeUnit_t      duration     = initTimeSeconds(duration_seconds);
    const angularUnit_t   inclination  = initAngleRadians(inclination_radians);
    const lengthUnit_t    distance     =  initLengthMpc(distance_mpc);
    const angularUnit_t   reference_orbital_phase = 
        initAngleRadians(reference_orbital_phase_in);
            
    // Setup companion structures:
    
    spin_t spin_1 = 
    {
        .x = spin_1_in[0],
        .y = spin_1_in[1],
        .z = spin_1_in[2]
    };
    spin_t spin_2 = 
    {
        .x = spin_2_in[0],
        .y = spin_2_in[1],
        .z = spin_2_in[2]
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
    
    const int32_t num_samples = 
        (int32_t)floor(sample_rate.hertz*duration.seconds);
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
    
    waveform_axes_s waveform_axes_td = 
        generateInspiral(
            system_properties,
            temporal_properties,
            num_waveforms,
            approximant
        );
    
    float2_t *strain = NULL;
    cudaToHost(
        (void**)&waveform_axes_td.strain.values[waveform_axes_td.strain.num_samples - num_samples - 1], 
        sizeof(float2_t),
        num_samples,
        (void**) &strain
    );
    cudaFree(waveform_axes_td.strain.values);
    cudaFree(waveform_axes_td.time.values);
        
    if (waveform_axes_td.strain.num_samples < num_samples) 
    {    
        fprintf(
            stderr, 
            "Warning! Cuphenom not generating waveforms of desired num_samples."
            "\n"
        );
    }
            
    float *test_array = malloc(sizeof(float) * (size_t)num_samples * 2);
    for (int32_t index = 0; index < num_samples; index++)
    {
        test_array[2*index + 0] = strain[index].x;
        test_array[2*index + 1] = strain[index].y;
    }
    
    free(strain);
        
    return test_array;
}

#endif