#ifndef PHENOM_STRUCTURES_H
#define PHENOM_STRUCTURES_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <cufft.h>
#include <cuda_fp16.h>

#include "phenomd_data.h"

// ~~~~~~~~~~ General structures ~~~~~~~~~~:

typedef struct 
{
    cuFloatComplex plus;
    cuFloatComplex cross;
} complex_strain_element_c;

typedef struct 
{
    complex float plus;
    complex float cross;
} complex_strain_element_t;

typedef struct
{
    frequencyUnit_t *values;
    frequencyUnit_t  interval;
    int32_t          num_samples;
} frequency_array_s;

typedef struct
{
    complex_strain_element_c *values;
    int32_t                   num_samples;
} complex_strain_array_s;

typedef struct
{
    frequency_array_s      frequency;
    complex_strain_array_s strain;
} complex_waveform_axes_s;

// Vector of dimensionless spins {x,y,z}:
typedef struct 
{ 
    float x;
    float y;
    float z;
} spin_t;

typedef struct 
{
    massUnit_t mass;           // <-- Mass of companion (massUnit_t).     
    spin_t spin;               // <-- Vector of dimensionless spins {x,y,z}.
    float quadrapole_moment;  
    float lambda;
} companion_s;

typedef struct 
{
    // ~~~~ Binary Companions ~~~~~ //
    companion_s companion[2];
    
    // ~~~~ Mass properties ~~~~~ //
    massUnit_t total_mass;
    massUnit_t reduced_mass;
    float      symmetric_mass_ratio;
    
    // ~~~~ Distance properties ~~~~~ //
    float       redshift;
    
    // Distance of source (lengthUnit_t)@
    lengthUnit_t distance;             

    // ~~~~ Orbital properties ~~~~~ //
    
    // Reference orbital phase (angularUnit_t):
    angularUnit_t reference_orbital_phase; 
    
    // Longitude of ascending nodes, degenerate with the polarization angle:
    float ascending_node_longitude;

    // Inclination of source (angularUnit_t):
    angularUnit_t inclination;              
    
    // Eccentrocity at reference epoch:
    float eccentricity;             

    // Mean anomaly of periastron:
    float mean_periastron_anomaly;  
} system_properties_s;

// Useful powers in GW waveforms: 1/6, 1/3, 2/3, 4/3, 5/3, 2, 7/3, 8/3, -1, 
// -1/6, -7/6, -1/3, -2/3, -5/3 calculated using only one invocation of 'pow', 
// the rest are just multiplications and divisions:
typedef struct {
    float one;
    float third;
    float two_thirds;
    float four_thirds;
    float five_thirds;
    float two;
    float seven_thirds;
    float eight_thirds;
    float inv;
    float m_seven_sixths;
    float m_third;
    float m_two_thirds;
    float m_five_thirds;
} useful_powers_s;

// ~~~~~~~~~~ Amplitude structures ~~~~~~~~~~:

typedef struct {
   float eta;         // symmetric mass-ratio
   float etaInv;      // 1/eta
   float chi12;       // chi1*chi1;
   float chi22;       // chi2*chi2;
   float eta2;        // eta*eta;
   float eta3;        // eta*eta*eta;
   float Seta;        // sqrt(1.0 - 4.0*eta);
   float SetaPlus1;   // (1.0 + Seta);
   float chi1, chi2;  // dimensionless aligned spins, convention m1 >= m2.
   float q;           // asymmetric mass-ratio (q>=1)
   float chi;         // PN reduced spin parameter
   float ringdown_frequency;         // ringdown frequency
   float damping_time;         // imaginary part of the ringdown frequency (damping time)
  
   // Frequency at which the mrerger-ringdown amplitude is maximum:
   float inspiral_merger_peak_frequency;    
  
   // Phenomenological inspiral amplitude coefficients:
   float inspiral[NUM_AMPLITUDE_INSPIRAL_COEFFICIENTS];
  
   // Phenomenological intermediate amplitude coefficients:
   float intermediate[NUM_AMPLITUDE_INTERMEDIATE_TERMS];
  
   // Phenomenological merger-ringdown amplitude coefficients:
   float merger_ringdown[NUM_AMPLITUDE_MERGER_RINGDOWN_COEFFICIENTS];
} amplitude_coefficients_s;

/**
 * used to cache the recurring (frequency-independent) prefactors of AmpInsAnsatz. Must be inited with a call to
 * init_amp_ins_prefactors(&prefactors, p);
 */
typedef struct {
     float two_thirds;
     float one;
     float four_thirds;
     float five_thirds;
     float two;
     float seven_thirds;
     float eight_thirds;
     float three;
     float amp0;
} amplitude_inspiral_prefactors_s;

// ~~~~~~~~~~ Phase structures ~~~~~~~~~~:

typedef struct {
   float eta;                // symmetric mass-ratio
   float etaInv;             // 1/eta
   float eta2;               // eta*eta
   float Seta;               // sqrt(1.0 - 4.0*eta);
   float chi1, chi2;         // dimensionless aligned spins, convention m1 >= m2.
   float q;                  // asymmetric mass-ratio (q>=1)
   float chi;                // PN reduced spin parameter
   float ringdown_frequency; // ringdown frequency
   float damping_time;       // imaginary part of the ringdown frequency (damping time)
  
   // Phenomenological inspiral phase coefficients
   float inspiral[NUM_PHASE_INSPIRAL_COEFFICIENTS];
    
   // Phenomenological intermediate phase coefficients
   float intermediate[NUM_PHASE_INTERMEDIATE_COEFFICIENTS];
  
   // Phenomenological merger-ringdown phase coefficients
   float merger_ringdown[NUM_PHASE_MERGER_RINGDOWN_COEFFICIENTS];
  
   // C1 phase connection coefficients
   float C1Int;
   float C2Int;
   float C1MRD;
   float C2MRD;
  
   // Transition frequencies:
   float intermediate_start_frequency;  
   float merger_ringdown_start_frequency;
} phase_coefficients_s;

typedef struct {
    float initial_phasing;
    float third;
    float third_with_logv;
    float two_thirds;
    float one;
    float four_thirds;
    float five_thirds;
    float two;
    float logv;
    float minus_third;
    float minus_two_thirds;
    float minus_one;
    float minus_four_thirds;
    float minus_five_thirds;
} phase_inspiral_prefactors_s;

// Structure for passing around PN phasing coefficients.
// For use with the TaylorF2 waveform:
#define PN_PHASING_SERIES_MAX_ORDER 15
typedef struct {
    float v      [PN_PHASING_SERIES_MAX_ORDER+1];
    float vlogv  [PN_PHASING_SERIES_MAX_ORDER+1];
    float vlogvsq[PN_PHASING_SERIES_MAX_ORDER+1];
} pn_phasing_series_s;

#endif



