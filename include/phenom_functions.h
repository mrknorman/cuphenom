#ifndef PHENOM_FUNCTIONS_H
#define PHENOM_FUNCTIONS_H

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include <cuda_maths.h>

float calculateSpinNorm(
    const spin_t spin
    ) {
    
    return spin.x*spin.x + spin.y*spin.y + spin.z*spin.z;
}

complex_waveform_axes_s generatePhenomD(
    const frequencyUnit_t  starting_frequency,
    const frequencyUnit_t  ending_frequency,
    const frequencyUnit_t  frequency_interval,
    const int32_t          num_strain_axis_samples
    );

int32_t sumPhenomDFrequencies(
          complex_waveform_axes_s          waveform_axes,
    const float                            inclination,
    const float                            total_mass_seconds,
    const amplitude_coefficients_s         amplitude_coefficients,
    const amplitude_inspiral_prefactors_s  amplitude_prefactors,
    const phase_coefficients_s             phase_coefficients, 
    const phase_inspiral_prefactors_s      phase_prefactors, 
    const int32_t                          offset,
    const float                            phase_shift,
    const float                            amp0,
    const float                            reference_mass_frequency,
    const float                            phi_precalc
    );

int32_t applyPolarization(
    const complex_waveform_axes_s waveform_axes,
    const float                   polarization
    );

int32_t taperWaveform(
    const complex_waveform_axes_s waveform_axes,
    const float                   starting_frequency,
    const float                   minimum_frequency,
    const float                   frequency_interval
    );

static void checkSystemParameters(
    const system_properties_s system_properties
    ) {
    
    // Unpack companion structs for readability:
    const companion_s companion_1 = system_properties.companion[0];
    const companion_s companion_2 = system_properties.companion[1];
    
    const massUnit_t   total_mass = system_properties.total_mass;
    const lengthUnit_t distance   = system_properties.distance;
    
    if (companion_1.mass.kilograms < 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! Mass 1 (%f) must be positive.\n",
            __func__,
            companion_1.mass.kilograms
        );   
    }
    if (companion_2.mass.kilograms < 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! Mass 2 (%f) must be positive.\n",
            __func__,
            companion_2.mass.kilograms
        );   
    }
    if (distance.meters <= 0.0f)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! Distance (%f) must be non-negative.\n",
            __func__,
            distance.meters
        );   
    }
    if (companion_1.mass.kilograms < initMassSolarMass(0.09f).kilograms)
    {
        fprintf(
            stderr, 
            "%s:\n Warning! Small value of m1 = %e (kg) = %e (Msun) requested.\n"
            "Perhaps you have a unit conversion error?\n", 
            __func__, 
            companion_1.mass.kilograms, 
            companion_1.mass.msun
        );
    }
    if (companion_2.mass.kilograms < initMassSolarMass(0.09f).kilograms)
    {
         fprintf(stderr,
            "%s:\n Small value of m2 = %e (kg) = %e (Msun) requested.\n"
            "Perhaps you have a unit conversion error?\n", 
            __func__, 
            companion_2.mass.kilograms, 
            companion_2.mass.msun
        );
    }
    if (total_mass.kilograms > initMassSolarMass(1000.0f).kilograms)
    {
         fprintf(
             stderr,
            "%s:\n Warning! Large value of total mass m1+m2 = %e (kg) = %e "
            "(Msun) requested.\nSignal not likely to be in band of ground-based"
            " detectors.\n", 
            __func__, 
            total_mass.kilograms, 
            total_mass.msun
        );
    }
    if (calculateSpinNorm(companion_1.spin) > 1.000001f)
    {
         fprintf(
             stderr, 
             "Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested.\nAre you "
             "sure you want to violate the Kerr bound?\n", 
             __func__, 
             companion_1.spin.x, companion_1.spin.y, companion_1.spin.z
        );
    }
    if(calculateSpinNorm(companion_2.spin) > 1.000001f)
    {
        fprintf(
            stderr, 
            "Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested.\nAre you "
            "sure you want to violate the Kerr bound?\n", 
            __func__, 
            companion_2.spin.x, companion_2.spin.y, companion_2.spin.z
        );
    }

    const float mass_ratio = 
        companion_1.mass.kilograms/companion_2.mass.kilograms;
    const float max_mass_ratio = 1000.0f;
    
    if (mass_ratio > max_mass_ratio)
    {    
        fprintf(
            stderr,
            "%s:\n"
            "Warning! Mass ratio (%f) is larger than maximum max ratio allowed "
            "for model (%f). \n",
            __func__,
            mass_ratio,
            max_mass_ratio
        );
    }
    if (
        (((companion_1.spin.z >  1.0f) || (companion_1.spin.z < -1.0f)) 
        || (companion_2.spin.z >  1.0f)) || (companion_2.spin.z < -1.0f)
    ) {
        fprintf(
           stderr,
           "%s:\n"
           "Warning! Spins outside the range [-1,1] are not supported. \n",
            __func__
        );
    }
}

system_properties_s initBinarySystem(
    companion_s   companion_a,
    companion_s   companion_b,
    lengthUnit_t  distance,
    float         redshift,
    angularUnit_t inclination,
    angularUnit_t reference_orbital_phase,
    float         ascending_node_longitude,
    float         eccentricity, 
    float         mean_periastron_anomaly
    ) {
    
    // Initilise structure to hold system_properties information:
    system_properties_s system_properties;
    
    // Set system_properties redshift:
    system_properties.redshift = redshift;

    // Apply redshift correction to dimensionful source-frame quantities:
    companion_a.mass = scaleMass(companion_a.mass, 1.0f + redshift);
    companion_b.mass = scaleMass(companion_b.mass, 1.0f + redshift);
    
    // Change from comoving (transverse) distance to luminosity distance:
    system_properties.distance = scaleLength(distance, 1.0f + redshift);  
    
    // Set companion one as higher mass input:
    if (companion_a.mass.kilograms > companion_b.mass.kilograms)
    {
        system_properties.companion[0] = companion_a;
        system_properties.companion[1] = companion_b;
    }
    else
    {
        system_properties.companion[1] = companion_a;
        system_properties.companion[0] = companion_b;
    }    
    
    // Calculate total mass of the system_properties:
    system_properties.total_mass = 
        addMasses(
            system_properties.companion[0].mass, 
            system_properties.companion[1].mass
        );
    
    // Calculate the reduced mass of the system_properties:
    system_properties.reduced_mass = 
        divideMasses(
            multiplyMasses(
                system_properties.companion[0].mass, 
                system_properties.companion[1].mass
            ), 
            system_properties.total_mass
        );
    
    // Calculate the symmetric mass ratio:
    system_properties.symmetric_mass_ratio = 
        system_properties.reduced_mass.msun / 
        system_properties.total_mass.msun;
    
    // Assign orbital properties:
    system_properties.inclination              = inclination;
    system_properties.reference_orbital_phase  = reference_orbital_phase;       
    system_properties.ascending_node_longitude = ascending_node_longitude;        
    system_properties.eccentricity             = eccentricity;
    system_properties.mean_periastron_anomaly  = mean_periastron_anomaly;    
    
    // General sanity check the system_properties parameters. This will only 
    // give warnings:
    checkSystemParameters(system_properties);
    
    return system_properties;
}

typedef struct {
    // Sampling interval (timeUnit_t):
    timeUnit_t        sampling_interval;  
    
    // Starting GW frequency (frequencyUnit_t):
    frequencyUnit_t   starting_frequency;  
    
    // Ending GW Frequency (frequencyUnit_t):
    frequencyUnit_t   ending_frequency;  
    
    // Reference GW frequency (frequencyUnit_t):
    frequencyUnit_t   reference_frequency; 
    
    // Extra time to include for all waveforms to take care of situations where 
    // the frequency is close to merger (and is sweeping rapidly) this is a few 
    // cycles at the low frequency:
    timeUnit_t extra_time;
    
    // Upper bound on the chirp time starting at starting_frequency:
    timeUnit_t chirp_time_upper_bound;
     
    // Upper bound on the plunge and merger time:
    timeUnit_t merge_time_upper_bound;
    
    // Upper bound on the ringdown time:
    timeUnit_t ringdown_time_upper_bound;
    
    // Upper bound on the total time:
    timeUnit_t total_time_upper_bound;
     
} temporal_properties_s;

void performTimeShiftHost_old(
          complex float   *h_plus_frequency, 
          complex float   *h_cross_frequency, 
    const float            temporal_properties,
    const int32_t          num_waveform_samples,
    const frequencyUnit_t  frequency_interval
);

int32_t performTimeShift(
          complex_waveform_axes_s waveform_axes,
    const float                   num_samples_time_shift
    );

inline float TaylorT2Timing_0PNCoeff(
    const massUnit_t total_mass,
    const float      sym_mass_ratio
    ) {
    return -5.0f*total_mass.seconds/(256.0f*sym_mass_ratio);
}

inline float TaylorT2Timing_2PNCoeff(
    const float sym_mass_ratio
    ) {
    return 7.43f/2.52f + 11.0f/3.0f * sym_mass_ratio;
}

inline float TaylorT2Timing_4PNCoeff(
    const float sym_mass_ratio
    ) {
    return 30.58673f/5.08032f + 54.29f/5.04f*sym_mass_ratio 
         + 61.7f/7.2f*sym_mass_ratio*sym_mass_ratio;
}

inline float TaylorT3Frequency_0PNCoeff(
    const massUnit_t mass
    ) {    
    return 1.0f / (8.0f*(float)M_PI*mass.seconds);
}

float InspiralFinalBlackHoleSpinBound(
    const system_properties_s system_properties
    ) {
    
    // Lower bound on the final plunge, merger, and ringdown time here the
    // final black hole spin is overestimated by using the formula in Tichy and
    // Marronetti, Physical Review D 78 081501 (2008), Eq. (1) and Table 1, for
    // equal mass black holes, or the larger of the two spins (which covers the
    // extreme mass case).

    // Function constants:
    const float maximum_black_hole_spin = 0.998f;
    
    // Unpack companion structs for readability:
    const spin_t spin_1 = system_properties.companion[0].spin;
    const spin_t spin_2 = system_properties.companion[1].spin;
    
    float final_spin_upper_bound = 0.686f + 0.15f*(spin_1.z + spin_2.z);
    final_spin_upper_bound = 
       (final_spin_upper_bound < fabsf(spin_1.z))*fabsf(spin_1.z) 
     + (final_spin_upper_bound > fabsf(spin_1.z))*final_spin_upper_bound;
    final_spin_upper_bound = 
       (final_spin_upper_bound < fabsf(spin_2.z))*fabsf(spin_2.z) 
     + (final_spin_upper_bound > fabsf(spin_2.z))*final_spin_upper_bound;

    // It is possible that |S1z| or |S2z| >= 1, but s must be less than 1
    // (0th law of thermodynamics) so provide a maximum value for s:
    final_spin_upper_bound = 
        (final_spin_upper_bound > maximum_black_hole_spin)
            * maximum_black_hole_spin
      + (final_spin_upper_bound < maximum_black_hole_spin)
            * final_spin_upper_bound;

     return final_spin_upper_bound;
}

temporal_properties_s fixReferenceFrequency(
          temporal_properties_s temporal_properties,
    const Approximant           approximant
) {
     if (temporal_properties.reference_frequency.hertz == 0.0f)
     {
        switch (approximant) 
        {
            case IMRPhenomXPHM:
                temporal_properties.reference_frequency = 
                    temporal_properties.starting_frequency;
            default:
                break;
        }
    }
    return temporal_properties;
}

timeUnit_t InspiralChirpTimeBound(
    const frequencyUnit_t        starting_frequency, 
    const system_properties_s    system_properties
) {
    
    // Unpack companion structs for readability:
    const spin_t spin_1 = system_properties.companion[0].spin;
    const spin_t spin_2 = system_properties.companion[1].spin;
    
    // Unpack properties for readability:
    const massUnit_t total_mass           = system_properties.total_mass;
    const float      symmetric_mass_ratio = 
        system_properties.symmetric_mass_ratio;
    
    // over-estimate of chi
    const float chi = fabsf(
            ((fabsf(spin_1.z) >  fabsf(spin_2.z))*spin_1.z)
         +  ((fabsf(spin_1.z) <= fabsf(spin_2.z))*spin_2.z)
         );
     
    const float c0 = 
        fabsf(
            TaylorT2Timing_0PNCoeff(
                total_mass, 
                symmetric_mass_ratio
            )
        );
    const float c2 = 
        TaylorT2Timing_2PNCoeff(symmetric_mass_ratio);
    
    // The 1.5pN spin term is in TaylorT2 is 8*beta/5
    // where beta = (113/12 + (25/4)(m2/m1))*(s1*m1^2/M^2) + 2 <-> 1
    // [Cutler & Flanagan, Physical Review D 49, 2658 (1994), Eq. (3.21)]
    // which can be written as (113/12)*chi - (19/6)(s1 + s2)
    // and we drop the negative contribution:
    const float c3 = (226.0f/15.0f) * chi;
     
    // There is also a 1.5PN term with eta, but it is negative so do not 
    // include it.
    const float c4 = TaylorT2Timing_4PNCoeff(symmetric_mass_ratio);
    const float v = 
        cbrtf((float)M_PI*G_SI*total_mass.kilograms*starting_frequency.hertz)/C_SI;
     
    return initTimeSeconds(
        c0 * powf(v, -8.0f) * (1.0f + (c2 + (c3 + c4 * v) * v) * v * v)
    );
}

inline timeUnit_t InspiralMergeTimeBound(
    const system_properties_s system_properties
) {
    
    // Unpack properties for readability:
    const massUnit_t total_mass = system_properties.total_mass;
    
    return initTimeSeconds(2.0f*(float)M_PI*((9.0f*total_mass.meters)/(C_SI/3.0f)));
}

timeUnit_t InspiralRingdownTimeBound(
    const system_properties_s system_properties
    ) {
    
    // Waveform generators only go up to 10:
    const float nefolds = 11.0f; 
    
    // Unpack properties for readability:
    const massUnit_t total_mass = system_properties.total_mass;

    // Upper bound on the final black hole spin:
    const float final_spin_upper_bound = 
        InspiralFinalBlackHoleSpinBound(system_properties);

    // These values come from Table VIII of Berti, Cardoso, and Will with n=0, 
    // m=2 :
    const float f[] = {1.5251f, -1.1568f,  0.1292f}; 
    const float q[] = {0.7000f,  1.4187f, -0.4990f}; 

    const float omega = 
          (f[0] + f[1]*powf(1.0f - final_spin_upper_bound, f[2]))
        / total_mass.seconds;
    const float Q = q[0] + q[1] * powf(1.0f - final_spin_upper_bound, q[2]);
    
    // See Eq. (2.1) of Berti, Cardoso, and Will:
    const float tau = 2.0f * Q / omega; 

    return initTimeSeconds(nefolds * tau);
}

frequencyUnit_t InspiralChirpStartFrequencyBound(
    const timeUnit_t          duration, 
    const system_properties_s system_properties
    ) {
    
    // Unpack properties for readability:
    const massUnit_t total_mass = 
        system_properties.total_mass;
    const float symmetric_mass_ratio = 
        system_properties.symmetric_mass_ratio;
     
    const float c0 = TaylorT3Frequency_0PNCoeff(total_mass);
    return initFrequencyHertz(
            c0*powf( 
                   5.0f * total_mass.seconds 
                / (symmetric_mass_ratio * duration.seconds), 
                  3.0f / 8.0f
            )
        );
}

static void checkFreqeuncyParameters(
    const timeUnit_t      sampling_interval,
    const frequencyUnit_t starting_frequency
    ) {
    
    if (sampling_interval.seconds > 1.0f)
    {
        fprintf(
            stderr,
            "Warning %s \n Large value of sampling_interval = %e (s) requested."
            "\nPerhaps sample rate and time step size were swapped?\n", 
            __func__, sampling_interval.seconds
        );
    }
    if (sampling_interval.seconds < 1.0f/16385.0f)
    {
        fprintf(
            stderr,
            "Warning %s \n Small value of sampling_interval = %e (s) requested."
            "\nCheck for errors, this could create very large time series.\n", 
            __func__, sampling_interval.seconds
        );
    }
    if(starting_frequency.hertz < 1.0f)
    {
        fprintf(
            stderr, 
            "Warning - %s: Small value of starting_frequency = %e (hz) "
            "requested.\n Check for errors, this could create a very long"
            "waveform.\n", 
            __func__, starting_frequency.hertz
        );
    }
    if(starting_frequency.hertz > 40.000001f)
    {
         fprintf(
             stderr, 
             "Warning - %s: Large value of starting_frequency = %e (hz)" 
             "requested.\n Check for errors, the signal will start in band.\n", 
             __func__, starting_frequency.hertz
        );
    }
}


inline frequencyUnit_t calculateKerrISCOFrequency(
    system_properties_s system_properties
    ) {
    
    return initFrequencyHertz(
        1.0f / (powf(9.0f, 1.5f)*(float)M_PI*system_properties.total_mass.seconds)
    );
}

temporal_properties_s initTemporalProperties(
          timeUnit_t          sampling_interval,   // <-- Sampling interval (timeUnit_t).
          frequencyUnit_t     starting_frequency,  // <-- Starting GW frequency (frequencyUnit_t).
          frequencyUnit_t     reference_frequency, // <-- Reference GW frequency (frequencyUnit_t).
    const system_properties_s system_properties,
    const Approximant         approximant
    ) {
    
    // General sanity check the temporal_properties parameters. This will only 
    // give warnings:
    checkFreqeuncyParameters(
       sampling_interval,
       starting_frequency
    );
    
    temporal_properties_s temporal_properties;
    
    temporal_properties.sampling_interval   = sampling_interval;
    temporal_properties.reference_frequency = reference_frequency;
    
    // Adjust the reference frequency for certain precessing approximants:
    // if that approximate interprets reference_frequency==0 to be 
    // starting_frequency, set reference_frequency=starting_frequency otherwise 
    // do nothing:
    temporal_properties = 
        fixReferenceFrequency(temporal_properties, approximant);
    
    // Calculate ending frequency:
    temporal_properties.ending_frequency = 
        initFrequencyHertz(0.5f/sampling_interval.seconds);

    // If the requested low frequency is below the lowest Kerr ISCO
    // frequency then change it to that frequency:
    const frequencyUnit_t kerr_isco_frequency = 
        calculateKerrISCOFrequency(system_properties);
    
    if (starting_frequency.hertz > kerr_isco_frequency.hertz)
         temporal_properties.starting_frequency = kerr_isco_frequency;
    else
        temporal_properties.starting_frequency = starting_frequency;
    
    // Calculate time boundaries:
    
    // Extra time to include for all waveforms to take care of situations where 
    // the frequency is close to merger (and is sweeping rapidly) this is a few 
    // cycles at the low frequency:
    temporal_properties.extra_time = 
        initTimeSeconds(
            EXTRA_CYCLES / temporal_properties.starting_frequency.hertz
        );

    // Upper bound on the chirp time starting at starting_frequency:
    temporal_properties.chirp_time_upper_bound =
        InspiralChirpTimeBound(
            temporal_properties.starting_frequency, 
            system_properties
        );

    // Upper bound on the plunge and merger time:
    temporal_properties.merge_time_upper_bound = 
        InspiralMergeTimeBound(system_properties);

    // Upper bound on the ringdown time:
    temporal_properties.ringdown_time_upper_bound = 
        InspiralRingdownTimeBound(system_properties);
    
    // Upper bound on the total time:
    temporal_properties.total_time_upper_bound = 
        addTimes(
            4,
            scaleTime(
                temporal_properties.chirp_time_upper_bound, 
                (1.0f + EXTRA_TIME_FRACTION)
            ),
            temporal_properties.merge_time_upper_bound,
            temporal_properties.ringdown_time_upper_bound,
            temporal_properties.extra_time
        );
        
    // Time domain approximant: condition by generating a waveform with a lower 
    // starting frequency and apply tapers in the region between that lower 
    // frequency and the requested frequency starting_frequency; here compute a 
    // new lower frequency:
    temporal_properties.starting_frequency = 
        InspiralChirpStartFrequencyBound(
            temporal_properties.total_time_upper_bound, 
            system_properties
        );
    
    return temporal_properties;
}

void castComplex64to32(
    const complex double  *in, 
          complex float   *out,  
    const         int32_t  num_elements
) {
    for(int32_t index = 0; index < num_elements; index++)
    {
        out[index] = (complex float) {creal(in[index]), cimag(in[index])};
    }
}

void cast32to64(float *in, double *out, int32_t n){
    for(int32_t i = 0; i< n; i++)
    {
        out[i] = ((double) in[i]);
    }
}

void performIRFFT(
          complex double                *hptilde,
          complex double                *hctilde,
                  double                *hplus,
                  double                *hcross,
                  temporal_properties_s  temporal_properties,
    const         frequencyUnit_t        frequency_interval,
    const         int32_t                num_waveform_samples
    ) {
    
    complex float  *float_array  = 
        (complex float*)malloc(
            sizeof(complex float)*(size_t)(num_waveform_samples*2)
        );    
    
    const int32_t num_frequency_bins = 
        (int32_t)round(num_waveform_samples/2) + 1; 
    castComplex64to32(
        hptilde,
        float_array, 
        num_frequency_bins
    );
    
    castComplex64to32(
        hctilde,
        &float_array[num_waveform_samples], 
        num_frequency_bins
    );

    complex float *ifft_g = NULL;
    cudaToDevice(
        (void*) float_array, 
        sizeof(complex float),
        num_waveform_samples*2,
        (void**) &ifft_g
    );    
    
    const float num_samples_time_shift = 
          roundf(temporal_properties.extra_time.seconds / 
          temporal_properties.sampling_interval.seconds) 
        * temporal_properties.sampling_interval.seconds;
    
    performTimeShiftHost_old(
          ifft_g, 
          &ifft_g[num_waveform_samples], 
          num_samples_time_shift,
          num_waveform_samples,
          frequency_interval
    );
    
    cudaIRfft(
        num_waveform_samples,
        2,
        (double) 
        ((float)num_waveform_samples * temporal_properties.sampling_interval.seconds) * 0.5, // I don't know where this factor of 0.5 comes from slightly concerning
        (cuFloatComplex*)ifft_g
    );
    
    float *ifft = NULL;
    cudaToHost(
        (void**)ifft_g, 
        sizeof(float),
        num_waveform_samples*3,
        (void**) &ifft
    );
    cast32to64(
        ifft, 
        hplus,
        num_waveform_samples
    );
    cast32to64(
        &ifft[2*num_waveform_samples], 
        hcross,
        num_waveform_samples
    );
    free(ifft);
    cudaFree(ifft_g);
}

void performIRFFT64(
    const complex double     *hptilde,
    const complex double     *hctilde,
                  double     *hplus,
                  double     *hcross,
                  timeUnit_t  sampling_interval,
    const         int32_t     num_waveform_samples
    ) 
    
    {
    complex double *float_array = 
        (complex double*)malloc(
            sizeof(complex double)*(size_t)(num_waveform_samples*2)
        );    
    
    for (int32_t index = 0; index < (num_waveform_samples/2 + 1); index++)
    {
        float_array[index]                        = hptilde[index];
        float_array[index + num_waveform_samples] = hctilde[index];
    }
    
    complex double *ifft_g = NULL;
    cudaToDevice(
        (void*) float_array, 
        sizeof(complex double),
        num_waveform_samples*2,
        (void**) &ifft_g
    );    
    
    cudaIRfft64(
        num_waveform_samples,
        2,
        (double)num_waveform_samples * sampling_interval.seconds,
        (cuDoubleComplex*)ifft_g
    );
    
    double *ifft = NULL;
    cudaToHost(
        (void**)ifft_g, 
        sizeof(double),
        num_waveform_samples*3,
        (void**) &ifft
    );
    
    for (int32_t index = 0; index < num_waveform_samples; index++)
    {
        hplus[index]  = ifft[index];
        hcross[index] = ifft[index + 2*num_waveform_samples];
    }
    
    free(ifft);
    cudaFree(ifft_g);
}

// The phasing function for TaylorF2 frequency-domain waveform.
// This function is tested in ../test/PNCoefficients.c for consistency
// with the energy and flux in this file.


// Taylor coefficients:

inline float TaylorF2Phasing_2PNCoeff(
    const float eta
    ) {
    return 5.0f*(74.3f/8.4f + 11.0f*eta)/9.0f;
}

inline float TaylorF2Phasing_3PNCoeff() {
    return -16.0f*(float)M_PI;
}
  
inline float TaylorF2Phasing_4PNCoeff(
    const float eta
    ) {
    return 5.0f*(3058.673f/7.056f + 5429.0f/7.0f*eta+617.0f*eta*eta)/72.0f;
}

inline float TaylorF2Phasing_5PNCoeff(
    const float eta
    ) {
    return 5.0f/9.0f*(772.9f/8.4f-13.0f*eta)*(float)M_PI;
}

inline float TaylorF2Phasing_5PNLogCoeff(
    const float eta
    ) {
    return 5.0f/3.0f*(772.9f/8.4f-13.0f*eta)*(float)M_PI;
}

inline float TaylorF2Phasing_6PNLogCoeff() {
    return -684.8f/2.1f;
}

inline float TaylorF2Phasing_6PNCoeff(
    const float eta
    ) {
   return 
         11583.231236531f/4.694215680f 
       - 640.0f/3.0f*PI_POWER_TWO 
       - 684.8f/2.1f*EULER_MASCHERONI 
       + eta*(-15737.765635f/3.048192f + 225.5f/1.2f*PI_POWER_TWO) 
       + eta*eta*76.055f/1.728f 
       - eta*eta*eta*127.825f/1.296f
       + TaylorF2Phasing_6PNLogCoeff()*logf(4.0f);
}

inline float TaylorF2Phasing_7PNCoeff(
    const float eta
    ) {
    return 
        (float)M_PI*(
              770.96675f/2.54016f 
            + 378.515f/1.512f*eta 
            - 740.45f/7.56f*eta*eta
        );
}

inline float TaylorF2Phasing_7PNSOCoeff(
    const float mByM
    ) {
    const float eta = mByM*(1.0f-mByM);
    
    return 
       mByM*(
           - 17097.8035f/4.8384f
           + eta*28764.25f/6.72f
           + eta*eta*47.35f/1.44f 
           + mByM*(
               - 7189.233785f/1.524096f
               + eta*458.555f/3.024f
               - eta*eta*534.5f/7.2f
            )
        );
}

// Tidal corrections to F2 phasing
// See arXiv:1101.1673
inline float TaylorF2Phasing_10PNTidalCoeff(
    const float mByM // < ratio of object mass to total mass
    ) {
    return (-288.0f + 264.0f*mByM)*mByM*mByM*mByM*mByM;
}

inline float TaylorF2Phasing_12PNTidalCoeff(
    const float mByM // Ratio of object mass to total mass
    ) {
    
    return 
        (- 15895.0f/28.0f
         + 4595.0f/28.0f*mByM 
         + 5715.0f/14.0f*mByM*mByM 
         - 325.0f/7.0f*mByM*mByM*mByM
        )*mByM*mByM*mByM*mByM;
}

inline float TaylorF2Phasing_13PNTidalCoeff(
    const float mByM /**< ratio of object mass to total mass */
    ) { 
    // literature: Agathos et al (arxiv 1503.0545) eq (5)
    // the coefficient mByM4 conversion & transformation (6.5PN, 7PN, 7.5PN):
    // mByM=mA/M: mA= mass star A, M is total mass (mA+mB)
    // Lambda (unitless) = lambda(m) / mA^5 
    // to call the function: 
    // Lambda * XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff 
    // lambda(m)*mByM^4/mA^5= lambda(m)*(mA/M)^4/(mA)^5= lambda/(M^4*mA) 
    // =lambda/(mByM*M^5) eq (5) 
    
    return mByM*mByM*mByM*mByM*24.0f*(12.0f - 11.0f*mByM)*(float)M_PI;
}
  

inline float InspiralTaylorF2Phasing_13PNTidalCoeff(
     const float mByM // ratio of object mass to total mass
     ) {
    //  literature: Agathos et al (arxiv 1503.0545) eq (5)
    // the coefficient mByM4 conversion & transformation (6.5PN, 7PN, 7.5PN):
    // mByM=mA/M: mA= mass star A, M is total mass (mA+mB)
    // Lambda (unitless) = lambda(m) / mA^5 
    // to call the function: 
    // Lambda * XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff 
    // lambda(m)*mByM^4/mA^5= lambda(m)*(mA/M)^4/(mA)^5= lambda/(M^4*mA) 
    // =lambda/(mByM*M^5) eq (5) 
    return mByM*mByM*mByM*mByM*24.0f*(12.0f - 11.0f*mByM)*(float)M_PI;
}

inline float TaylorF2Phasing_14PNTidalCoeff(
    const float mByM // ratio of object mass to total mass
    ) {
    //literature: Agathos et al (arxiv 1503.0545) eq (5)
    //caveat: these are incomplete terms
    //conversion see XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff above
    //--> completed by the terms given in equation (4) of :
    //Tatsuya Narikawa, Nami Uchikata, Takahiro Tanaka,
    //"Gravitational-wave constraints on the GWTC-2 events by measuring
    //the tidal deformability and the spin-induced quadrupole moment",
    //Phys. Rev. D 104, 084056 (2021), arXiv:2106.09193
    const float mByM3 = mByM*mByM*mByM;
    const float mByM4 = mByM3 * mByM;
    return 
        - mByM4*5.0f*(
              193986935.0f/571536.0f 
            - 14415613.0f/381024.0f*mByM 
            - 57859.0f/378.0f*mByM*mByM 
            - 209495.0f/1512.0f*mByM3 
            + 965.0f/54.0f*mByM4 
            - 4.00f*mByM4*mByM
        );
}

inline float TaylorF2Phasing_15PNTidalCoeff(
    const float mByM // Ratio of object mass to total mass
    ) {
    //literature: Agathos et al (arxiv 1503.0545) eq (5)
    //conversion see XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff above 
    //--> corrected by the terms given in equation (4) of :
    //Tatsuya Narikawa, Nami Uchikata, Takahiro Tanaka,
    //"Gravitational-wave constraints on the GWTC-2 events by measuring
    //the tidal deformability and the spin-induced quadrupole moment",
    //Phys. Rev. D 104, 084056 (2021), arXiv:2106.09193
    
    const float mByM2 = mByM*mByM;
    const float mByM3 = mByM2*mByM;
    const float mByM4 = mByM3*mByM;
    return 
        mByM4*1.0f/28.0f*(float)M_PI*(
              27719.0f 
            - 22415.0f*mByM 
            + 7598.0f*mByM2 
            - 10520.0f*mByM3
        );
}

pn_phasing_series_s PNPhasing_F2(
    const float   m1, // Masns of body 1, in Msol
    const float   m2, // Mass of body 2, in Msol
    const float   chi1L, // Component of dimensionless spin 1 along Lhat
    const float   chi2L, // Component of dimensionless spin 2 along Lhat
    const float   lambda1,
    const float   lambda2,
    const int32_t tidal_pn_order
    ) {
    
    const float mtot = m1 + m2;
    const float eta  = m1*m2/mtot/mtot;
    const float m1M  = m1/mtot;
    const float m2M  = m2/mtot;

    const float pfaN = 3.0f/(128.0f * eta);

    pn_phasing_series_s pfa;
    
    memset(pfa.v, 0, sizeof(double) * (PN_PHASING_SERIES_MAX_ORDER + 1));
    
    pfa.v[0]     = 1.0f;
    pfa.v[1]     = 0.0f;
    pfa.v[2]     = TaylorF2Phasing_2PNCoeff(eta);
    pfa.v[3]     = TaylorF2Phasing_3PNCoeff();
    pfa.v[4]     = TaylorF2Phasing_4PNCoeff(eta);
    pfa.v[5]     = TaylorF2Phasing_5PNCoeff(eta);
    pfa.vlogv[5] = TaylorF2Phasing_5PNLogCoeff(eta);
    pfa.v[6]     = TaylorF2Phasing_6PNCoeff(eta);
    pfa.vlogv[6] = TaylorF2Phasing_6PNLogCoeff();
    pfa.v[7]     = TaylorF2Phasing_7PNCoeff(eta)
                 + TaylorF2Phasing_7PNSOCoeff(m1M)*chi1L
                 + TaylorF2Phasing_7PNSOCoeff(m2M)*chi2L;
                 
    switch(tidal_pn_order)
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_75PN:
            pfa.v[15] = lambda1*TaylorF2Phasing_15PNTidalCoeff(m1M) 
                      + lambda2*TaylorF2Phasing_15PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_7PN:
            pfa.v[14] = lambda1*TaylorF2Phasing_14PNTidalCoeff(m1M) 
                      + lambda2*TaylorF2Phasing_14PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_65PN:
            pfa.v[13] = lambda1*TaylorF2Phasing_13PNTidalCoeff(m1M) 
                      + lambda2*TaylorF2Phasing_13PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            pfa.v[12] = lambda1*TaylorF2Phasing_12PNTidalCoeff(m1M) 
                      + lambda2*TaylorF2Phasing_12PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pfa.v[10] = lambda1*TaylorF2Phasing_10PNTidalCoeff(m1M) 
                      + lambda2*TaylorF2Phasing_10PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
        fprintf(
            stderr, 
            "%s:\nWarning! Invalid tidal PN order (%i)\n.",
            __func__, tidal_pn_order
        );    
        break;
    }
    
    // At the very end, multiply everything in the series by pfaN:
    for(int32_t index = 0; index < PN_PHASING_SERIES_MAX_ORDER + 1; index++)
    {
        pfa.v      [index] *= pfaN;
        pfa.vlogv  [index] *= pfaN;
        pfa.vlogvsq[index] *= pfaN;
    }
    
    return pfa;
}

// brief Returns structure containing TaylorF2 phasing coefficients for given
// physical parameters.
pn_phasing_series_s initTaylorF2AlignedPhasingSeries(
    const system_properties_s system_properties,
    const int32_t             tidal_pn_order
    ) {
    
    const companion_s companion_1 = system_properties.companion[0];
    const companion_s companion_2 = system_properties.companion[1];
    
    pn_phasing_series_s phasing_series =
        PNPhasing_F2(
            companion_1.mass.msun, 
            companion_2.mass.msun, 
            companion_1.spin.z, 
            companion_2.spin.z, 
            companion_1.lambda,
            companion_2.lambda,
            tidal_pn_order
        );

    return phasing_series;
}


//XLALSimInspiralTaylorF2AlignedPhasing
  

#endif