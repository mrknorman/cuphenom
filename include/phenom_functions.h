#ifndef PHENOM_FUNCTIONS_H
#define PHENOM_FUNCTIONS_H

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include <cuda_maths.h>

// Vector of dimensionless spins {x,y,z}:
typedef struct 
{ 
    double x;
    double y;
    double z;
} spin_t;

float64_t calculateSpinNorm(
    const spin_t spin
    ) {
    
    return spin.x*spin.x + spin.y*spin.y + spin.z*spin.z;
}

typedef struct {
    massUnit_t mass;               // <-- Mass of companion (massUnit_t).     
    spin_t spin;               // <-- Vector of dimensionless spins {x,y,z}.
    double quadrapole_moment;  
    double lambda;
} companion_s;

typedef struct {
    
    // ~~~~ Binary Companions ~~~~~ //
    companion_s companion[2];
    
    // ~~~~ Mass properties ~~~~~ //
    massUnit_t total_mass;
    massUnit_t reduced_mass;
    double symmetric_mass_ratio;
    
    // ~~~~ Distance properties ~~~~~ //
    double       redshift;
    
    // Distance of source (lengthUnit_t)@
    lengthUnit_t distance;             

    // ~~~~ Orbital properties ~~~~~ //
    
    // Reference orbital phase (angularUnit_t):
    angularUnit_t reference_orbital_phase; 
    
    // Longitude of ascending nodes, degenerate with the polarization angle:
    double ascending_node_longitude;

    // Inclination of source (angularUnit_t):
    angularUnit_t inclination;              
    
    // Eccentrocity at reference epoch:
    double eccentricity;             

    // Mean anomaly of periastron:
    double mean_periastron_anomaly;  
} system_properties_s;

static void checkSystemParameters(
    const system_properties_s system_properties
    ) {
    
    // Unpack companion structs for readability:
    const companion_s companion_1 = system_properties.companion[0];
    const companion_s companion_2 = system_properties.companion[1];
    
    const massUnit_t total_mass = system_properties.total_mass;
    
    if (companion_1.mass.kilograms < initMassSolarMass(0.09).kilograms)
    {
        fprintf(
            stderr, 
            "Warning %s \n Small value of m1 = %e (kg) = %e (Msun) requested.\n"
            " Perhaps you have a unit conversion error?\n", 
            __func__, companion_1.mass.kilograms, companion_1.mass.msun
        );
    }
    if (companion_2.mass.kilograms < initMassSolarMass(0.09).kilograms)
    {
         fprintf(stderr,
            "%s:\n Small value of m2 = %e (kg) = %e (Msun) requested.\n"
            "Perhaps you have a unit conversion error?\n", 
            __func__, companion_2.mass.kilograms, companion_2.mass.msun
        );
    }
    if (total_mass.kilograms > initMassSolarMass(1000).kilograms)
    {
         fprintf(
             stderr,
            "%s:\n Warning! Large value of total mass m1+m2 = %e (kg) = %e "
            "(Msun) requested.\nSignal not likely to be in band of ground-based"
            " detectors.\n", 
            __func__, total_mass.kilograms, total_mass.msun
        );
    }
    if (calculateSpinNorm(companion_1.spin) > 1.000001)
    {
         fprintf(
             stderr, 
             "Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested.\nAre you "
             "sure you want to violate the Kerr bound?\n", 
             __func__, 
             companion_1.spin.x, companion_1.spin.y, companion_1.spin.z
        );
    }
    if(calculateSpinNorm(companion_2.spin) > 1.000001)
    {
        fprintf(
            stderr, 
            "Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested.\nAre you "
            "sure you want to violate the Kerr bound?\n", 
            __func__, 
            companion_2.spin.x, companion_2.spin.y, companion_2.spin.z
        );
    }
}

system_properties_s initBinarySystem(
    companion_s   companion_a,
    companion_s   companion_b,
    lengthUnit_t  distance,
    double        redshift,
    angularUnit_t inclination,
    angularUnit_t reference_orbital_phase,
    double        ascending_node_longitude,
    double        eccentricity, 
    double        mean_periastron_anomaly
    ) {
    
    // Initilise structure to hold system_properties information:
    system_properties_s system_properties;
    
    // Set system_properties redshift:
    system_properties.redshift = redshift;

    // Apply redshift correction to dimensionful source-frame quantities:
    companion_a.mass = scaleMass(companion_a.mass, 1.0 + redshift);
    companion_b.mass = scaleMass(companion_b.mass, 1.0 + redshift);
    
    // Change from comoving (transverse) distance to luminosity distance:
    system_properties.distance = scaleLength(distance, 1.0 + redshift);  
    
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
        system_properties.reduced_mass.kilograms / 
        system_properties.total_mass.kilograms;
    
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

void performTimeShiftHost(
          complex float        *h_plus_frequency, 
          complex float        *h_cross_frequency, 
    const temporal_properties_s temporal_properties,
    const int32_t               num_waveform_samples,
    const frequencyUnit_t       frequency_interval
);

inline float64_t TaylorT2Timing_0PNCoeff(
    const massUnit_t total_mass,
    const float64_t  sym_mass_ratio
    ) {
    return -5.0*total_mass.seconds/(256.0*sym_mass_ratio);
}

inline float64_t TaylorT2Timing_2PNCoeff(
    const float64_t sym_mass_ratio
    ) {
    return 7.43/2.52 + 11./3. * sym_mass_ratio;
}

inline float64_t TaylorT2Timing_4PNCoeff(
    const float64_t sym_mass_ratio
    ) {
    return 30.58673/5.08032 + 54.29/5.04*sym_mass_ratio 
         + 61.7/7.2*sym_mass_ratio*sym_mass_ratio;
}

inline double TaylorT3Frequency_0PNCoeff(
    const massUnit_t mass
    ) {    
    return 1.0 / (8.0*M_PI*mass.seconds);
}

float64_t InspiralFinalBlackHoleSpinBound(
    const system_properties_s system_properties
    ) {
    
    // Lower bound on the final plunge, merger, and ringdown time here the
    // final black hole spin is overestimated by using the formula in Tichy and
    // Marronetti, Physical Review D 78 081501 (2008), Eq. (1) and Table 1, for
    // equal mass black holes, or the larger of the two spins (which covers the
    // extreme mass case).

    // Function constants:
    const float64_t maximum_black_hole_spin = 0.998;
    
    // Unpack companion structs for readability:
    const spin_t spin_1 = system_properties.companion[0].spin;
    const spin_t spin_2 = system_properties.companion[1].spin;
    
    float64_t final_spin_upper_bound = 0.686 + 0.15 * (spin_1.z + spin_2.z);
    final_spin_upper_bound = 
       (final_spin_upper_bound < fabs(spin_1.z))*fabs(spin_1.z) 
     + (final_spin_upper_bound > fabs(spin_1.z))*final_spin_upper_bound;
    final_spin_upper_bound = 
       (final_spin_upper_bound < fabs(spin_2.z))*fabs(spin_2.z) 
     + (final_spin_upper_bound > fabs(spin_2.z))*final_spin_upper_bound;

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
     if (temporal_properties.reference_frequency.hertz == 0)
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
    const double     symmetric_mass_ratio = 
        system_properties.symmetric_mass_ratio;
    
    // over-estimate of chi
    const float64_t chi = fabs(
            ((fabs(spin_1.z) >  fabs(spin_2.z))*spin_1.z)
         +  ((fabs(spin_1.z) <= fabs(spin_2.z))*spin_2.z)
         );
     
    const float64_t c0 = 
        fabs(
            TaylorT2Timing_0PNCoeff(
                total_mass, 
                symmetric_mass_ratio
            )
        );
    const float64_t c2 = 
        TaylorT2Timing_2PNCoeff(symmetric_mass_ratio);
    
    // The 1.5pN spin term is in TaylorT2 is 8*beta/5
    // where beta = (113/12 + (25/4)(m2/m1))*(s1*m1^2/M^2) + 2 <-> 1
    // [Cutler & Flanagan, Physical Review D 49, 2658 (1994), Eq. (3.21)]
    // which can be written as (113/12)*chi - (19/6)(s1 + s2)
    // and we drop the negative contribution:
    const float64_t c3 = (226.0/15.0) * chi;
     
    // There is also a 1.5PN term with eta, but it is negative so do not 
    // include it.
    const float64_t c4 = TaylorT2Timing_4PNCoeff(symmetric_mass_ratio);
    const float64_t v = 
        cbrt(M_PI*G_SI*total_mass.kilograms*starting_frequency.hertz)/C_SI;
     
    return initTimeSeconds(
        c0 * pow(v, -8) * (1.0 + (c2 + (c3 + c4 * v) * v) * v * v)
    );
}

inline timeUnit_t InspiralMergeTimeBound(
    const system_properties_s system_properties
) {
    
    // Unpack properties for readability:
    const massUnit_t total_mass = system_properties.total_mass;
    
    return initTimeSeconds(2.0*M_PI*((9.0*total_mass.meters)/(C_SI/3.0)));
}

timeUnit_t InspiralRingdownTimeBound(
    const system_properties_s system_properties
    ) {
    
    // Waveform generators only go up to 10:
    const float64_t nefolds = 11; 
    
    // Unpack properties for readability:
    const massUnit_t total_mass = system_properties.total_mass;

    // Upper bound on the final black hole spin:
    const float64_t final_spin_upper_bound = 
        InspiralFinalBlackHoleSpinBound(system_properties);

    // These values come from Table VIII of Berti, Cardoso, and Will with n=0, 
    // m=2 :
    const float64_t f[] = {1.5251, -1.1568,  0.1292}; 
    const float64_t q[] = {0.7000,  1.4187, -0.4990}; 

    const float64_t omega = 
          (f[0] + f[1]*pow(1.0 - final_spin_upper_bound, f[2]))
        / total_mass.seconds;
    const float64_t Q = q[0] + q[1] * pow(1.0 - final_spin_upper_bound, q[2]);
    
    // See Eq. (2.1) of Berti, Cardoso, and Will:
    const float64_t tau = 2.0 * Q / omega; 

    return initTimeSeconds(nefolds * tau);
}

frequencyUnit_t InspiralChirpStartFrequencyBound(
    const timeUnit_t          duration, 
    const system_properties_s system_properties
    ) {
    
    // Unpack properties for readability:
    const massUnit_t total_mass = 
        system_properties.total_mass;
    const double     symmetric_mass_ratio = 
        system_properties.symmetric_mass_ratio;
     
    double c0 = TaylorT3Frequency_0PNCoeff(total_mass);
    return initFrequencyHertz(
            c0*pow( 
                   5.0 * total_mass.seconds 
                / (symmetric_mass_ratio * duration.seconds), 
                  3.0 / 8.0
            )
        );
}

static void checkFreqeuncyParameters(
    const timeUnit_t      sampling_interval,
    const frequencyUnit_t starting_frequency
    ) {
    
    if (sampling_interval.seconds > 1.0)
    {
        fprintf(
            stderr,
            "Warning %s \n Large value of sampling_interval = %e (s) requested."
            "\nPerhaps sample rate and time step size were swapped?\n", 
            __func__, sampling_interval.seconds
        );
    }
    if (sampling_interval.seconds < 1.0/16385.0)
    {
        fprintf(
            stderr,
            "Warning %s \n Small value of sampling_interval = %e (s) requested."
            "\nCheck for errors, this could create very large time series.\n", 
            __func__, sampling_interval.seconds
        );
    }
    if(starting_frequency.hertz < 1.0)
    {
        fprintf(
            stderr, 
            "Warning - %s: Small value of starting_frequency = %e (hz) "
            "requested.\n Check for errors, this could create a very long"
            "waveform.\n", 
            __func__, starting_frequency.hertz
        );
    }
    if(starting_frequency.hertz > 40.000001)
    {
         fprintf(
             stderr, 
             "Warning - %s: Large value of starting_frequency = %e (hz)" 
             "requested.\n Check for errors, the signal will start in band.\n", 
             __func__, starting_frequency.hertz
        );
    }
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
        initFrequencyHertz(0.5/sampling_interval.seconds);

    // If the requested low frequency is below the lowest Kerr ISCO
    // frequency then change it to that frequency:
    const frequencyUnit_t kerr_isco_frequency = 
        initFrequencyHertz(
            1.0 / (pow(9.0, 1.5)*M_PI*system_properties.total_mass.seconds)
        );
    
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
                (1.0 + EXTRA_TIME_FRACTION)
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
    for(int32_t i = 0; i< n; i++){
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
    
    performTimeShiftHost(
          ifft_g, 
          &ifft_g[num_waveform_samples], 
          temporal_properties,
          num_waveform_samples,
          frequency_interval
    );
    
    cudaIRfft(
        num_waveform_samples,
        2,
        (double) (num_waveform_samples * temporal_properties.sampling_interval.seconds) * 0.5, // I don't know where this factor of 0.5 comes from slightly concerning
        (cuFloatComplex*)ifft_g
    );
    
    float *ifft = NULL;
    cudaToHost(
        ifft_g, 
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
        ifft_g, 
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

#endif