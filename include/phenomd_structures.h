#ifndef PHENOM_STRUCTURES_H
#define PHENOM_STRUCTURES_H

#include "phenomd_data.h"

// ~~~~~~~~~~ General structures ~~~~~~~~~~:

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
   double eta;         // symmetric mass-ratio
   double etaInv;      // 1/eta
   double chi12;       // chi1*chi1;
   double chi22;       // chi2*chi2;
   double eta2;        // eta*eta;
   double eta3;        // eta*eta*eta;
   double Seta;        // sqrt(1.0 - 4.0*eta);
   double SetaPlus1;   // (1.0 + Seta);
   double chi1, chi2;  // dimensionless aligned spins, convention m1 >= m2.
   double q;           // asymmetric mass-ratio (q>=1)
   double chi;         // PN reduced spin parameter
   double fRD;         // ringdown frequency
   double fDM;         // imaginary part of the ringdown frequency (damping time)
  
   double fmaxCalc;    // frequency at which the mrerger-ringdown amplitude is maximum
  
   // Phenomenological inspiral amplitude coefficients
   double inspiral[NUM_AMPLITUDE_INSPIRAL_COEFFICIENTS];
  
   // Phenomenological intermediate amplitude coefficients
   double intermediate[NUM_DELTA_TERMS];
  
   // Phenomenological merger-ringdown amplitude coefficients
   double merger_ringdown[NUM_AMPLITUDE_MERGER_RINGDOWN_COEFFICIENTS];
  
   // Coefficients for collocation method. Used in intermediate amplitude model
   double f1, f2, f3;
   double v1, v2, v3;
   double d1, d2;
  
   // Transition frequencies for amplitude
   // We don't *have* to store them, but it may be clearer.
   double fInsJoin;    // Ins = Inspiral
   double fMRDJoin;    // MRD = Merger-Ringdown
} IMRPhenomDAmplitudeCoefficients;


/**
 * used to cache the recurring (frequency-independent) prefactors of AmpInsAnsatz. Must be inited with a call to
 * init_amp_ins_prefactors(&prefactors, p);
 */
typedef struct {
     double two_thirds;
     double one;
     double four_thirds;
     double five_thirds;
     double two;
     double seven_thirds;
     double eight_thirds;
     double three;
     double amp0;
} AmpInsPrefactors;

// ~~~~~~~~~~ Phase structures ~~~~~~~~~~:

typedef struct {
   double eta;         // symmetric mass-ratio
   double etaInv;      // 1/eta
   double eta2;        // eta*eta
   double Seta;        // sqrt(1.0 - 4.0*eta);
   double chi1, chi2;  // dimensionless aligned spins, convention m1 >= m2.
   double q;           // asymmetric mass-ratio (q>=1)
   double chi;         // PN reduced spin parameter
   double fRD;         // ringdown frequency
   double fDM;         // imaginary part of the ringdown frequency (damping time)
  
   // Phenomenological inspiral phase coefficients
   double sigma1;
   double sigma2;
   double sigma3;
   double sigma4;
   double sigma5;
  
   // Phenomenological intermediate phase coefficients
   double beta1;
   double beta2;
   double beta3;
  
   // Phenomenological merger-ringdown phase coefficients
   double alpha1;
   double alpha2;
   double alpha3;
   double alpha4;
   double alpha5;
  
   // C1 phase connection coefficients
   double C1Int;
   double C2Int;
   double C1MRD;
   double C2MRD;
  
   // Transition frequencies for phase
   double fInsJoin;    // Ins = Inspiral
   double fMRDJoin;    // MRD = Merger-Ringdown
} IMRPhenomDPhaseCoefficients;

typedef struct {
    double initial_phasing;
    double third;
    double third_with_logv;
    double two_thirds;
    double one;
    double four_thirds;
    double five_thirds;
    double two;
    double logv;
    double minus_third;
    double minus_two_thirds;
    double minus_one;
    double minus_four_thirds;
    double minus_five_thirds;
} PhiInsPrefactors;


// Structure for passing around PN phasing coefficients.
// For use with the TaylorF2 waveform:
#define PN_PHASING_SERIES_MAX_ORDER 15
typedef struct {
    float v      [PN_PHASING_SERIES_MAX_ORDER+1];
    float vlogv  [PN_PHASING_SERIES_MAX_ORDER+1];
    float vlogvsq[PN_PHASING_SERIES_MAX_ORDER+1];
} pn_phasing_series_s;

#endif



