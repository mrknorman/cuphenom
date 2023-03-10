#include "phenom_d_data.h"
#include "phenom_d_structures.h"

inline float Square(
    const float number
) {
    return number*number;
}

inline float Cube(
    const float number)
{
    return number*number*number;
}

inline float Quart(float number)
{
    float pow2 = Square(number);
    return pow2 * pow2;
}

inline int32_t calcNextPow2(const float n)
{
   // Use pow here, not bit-wise shift, as the latter seems to run against an
   // upper cutoff long before SIZE_MAX, at least on some platforms:
   return (int32_t) powf(2.0f,ceilf(log2f(n)));
}

// See phenomd_data.h for the coefficient terms.
inline float calcCoefficient(
    const float *terms,
    const float  eta,
    const float  eta2, 
    const float  xi
    ) {
    
   return 
      terms[0] + terms[1]*eta
   + (terms[2] + terms[3]*eta + terms[ 4]*eta2
   + (terms[5] + terms[6]*eta + terms[ 7]*eta2)*xi
   + (terms[8] + terms[9]*eta + terms[10]*eta2)*xi*xi)*xi;
}

void calcCoefficients(
    const float  **terms,
    const size_t   num_terms,
    const float    eta,
    const float    eta2, 
    const float    xi,
          float   *coeffcients
    ) {
    
    for (size_t index = 0; index < num_terms; index++)
    {
        coeffcients[index] = 
            calcCoefficient(
                terms[index],
                eta, 
                eta2, 
                xi
            );
    }
}
 
/*************************** Amplitude functions ******************************/

//////////////////// Amplitude: Merger-Ringdown functions //////////////////////


// Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
// (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
// was not available when PhenomD was tuned.
inline float subtract3PNSS(
    const system_properties_s system_properties
    ) {
    
    const float mass_1_msun = system_properties.companion[0].mass.msun;
    const float mass_2_msun = system_properties.companion[1].mass.msun;
    
    const float m1M = mass_1_msun / system_properties.total_mass.msun;
    const float m2M = mass_2_msun / system_properties.total_mass.msun;
    
    const float chi_1 = system_properties.companion[0].spin.z;
    const float chi_2 = system_properties.companion[1].spin.z;
    
    const float symmetric_mass_ratio = system_properties.symmetric_mass_ratio;
   
    return 
        (326.75f/1.12f + 557.5f/1.8f*symmetric_mass_ratio)*(symmetric_mass_ratio*chi_1*chi_2)
        + (
            + (4703.5f/8.4f+2935.0f/6.0f*m1M-120.0f*m1M*m1M) 
            + (-4108.25f/6.72f-108.5f/1.2f*m1M+125.5f/3.6f*m1M*m1M)
        )*m1M*m1M*(chi_1*chi_1)
        + (
            + (4703.5f/8.4f+2935.0f/6.0f*m2M-120.0f*m2M*m2M) 
            + (-4108.25f/6.72f-108.5f/1.2f*m2M+125.5f/3.6f*m2M*m2M)
        )*m2M*m2M*(chi_2*chi_2);
}

float chiPN(
    const float Seta, 
    const float eta, 
    const float chi1, 
    const float chi2
) {
   // Convention m1 >= m2 and chi1 is the spin on m1
   // The 0.5 factor missing in the definitions of chi_s and chi_a is
   // recovered in the return expresion:
   const float chi_s = (chi1 + chi2);
   const float chi_a = (chi1 - chi2);
  
   return 0.5f * (chi_s * (1.0f - eta * 76.0f / 113.0f) + Seta * chi_a);
}

// Formula to predict the total radiated energy. Equation 3.7 and 3.8 
// arXiv:1508.07250.
// Input parameter s defined around Equation 3.7 and 3.8.
inline float _predictIrradiatedEnergy(
    const float eta, 
    const float s
    ) {
    
    const float eta_2 = eta * eta;
    const float eta_3 = eta_2 * eta;
  
     return
         (
             eta*
             (
                   0.055974469826360077f 
                 + 0.5809510763115132f*eta 
                 - 0.9606726679372312f*eta_2 
                 + 3.352411249771192f*eta_3
            )*(
                1.0f
                +(
                    - 0.0030302335878845507f 
                    - 2.0066110851351073f*eta 
                    + 7.7050567802399215f*eta_2
                )*s
            )
        )/(
            1.0f 
            + (
                - 0.6714403054720589f 
                - 1.4756929437702908f*eta 
                + 7.304676214885011f*eta_2
            )*s
        );
}

float predictIrradiatedEnergy(
    const float eta, 
    const float chi1, 
    const float chi2
    ) {
     // Convention m1 >= m2
     const float Seta = sqrtf(1.0f - 4.0f * eta);
     const float m1   = 0.5f * (1.0f + Seta);
     const float m2   = 0.5f * (1.0f - Seta);
     const float m1s   = m1 * m1;
     const float m2s   = m2 * m2;
     // arXiv:1508.07250
     const float s = (m1s * chi1 + m2s * chi2) / (m1s + m2s);
  
     return _predictIrradiatedEnergy(eta, s);
}

float calculateRingdownFrequency(
    const float   eta, 
    const float   chi1, 
    const float   chi2, 
    const float   final_spin,
    const double *QNM_data
    ) {

    if (final_spin > 1.0f)
    {
        fprintf(
            stderr,
            "%s:\n"
            "Warning! Final spin (%f) > 1.0 not supported!\n",
            __func__, final_spin
            );
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *i_ringdown_frequency = 
       gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
    gsl_spline_init(
       i_ringdown_frequency, QNMData_a, QNM_data, QNMData_length
    );

    const float return_value = 
       (float)gsl_spline_eval(i_ringdown_frequency, (double)final_spin, acc) 
       / (1.0f - predictIrradiatedEnergy(eta, chi1, chi2));

    gsl_spline_free(i_ringdown_frequency);
    gsl_interp_accel_free(acc);

    return return_value;
}

typedef enum tagNRTidal_version_type {
  NRTidal_V, // version NRTidal: based on https://arxiv.org/pdf/1706.02969.pdf
  NRTidalv2_V, // version NRTidalv2: https://arxiv.org/abs/1905.06011
  NRTidalv2NoAmpCorr_V, //<version NRTidalv2, without amplitude corrections
  NRTidalv2NSBH_V, // version NRTidalv2: https://arxiv.org/abs/1905.06011 with 
                   // amplitude corrections for NSBH 
                   // (used for SEOBNRv4ROM_NRTidalv2_NSBH)
  NoNRT_V          // special case for PhenomPv2 BBH baseline */
} NRTidal_version_type;

///////////////////// Amplitude: Intermediate functions ////////////////////////
 
// The Newtonian term in LAL is fine and we should use exactly the same (either 
// hardcoded or call). We just use the Mathematica expression for convenience.

// Inspiral amplitude plus rho phenom coefficents. rho coefficients computed
// in rho1_fun, rho2_fun, rho3_fun functions. Amplitude is a re-expansion. 
// See 1508.07253 and Equation 29, 30 and Appendix B arXiv:1508.07253 for 
// details:
inline float calculateInspiralAmplitudeAnsatz(
    const useful_powers_s  mass_frequency, 
    const amplitude_inspiral_prefactors_s prefactors
    ) {
    
    return 1.0f 
        + mass_frequency.two_thirds   * prefactors.two_thirds
        + mass_frequency.four_thirds  * prefactors.four_thirds
        + mass_frequency.five_thirds  * prefactors.five_thirds
        + mass_frequency.seven_thirds * prefactors.seven_thirds 
        + mass_frequency.eight_thirds * prefactors.eight_thirds
        + mass_frequency.one 
        * (   
            + prefactors.one 
            + mass_frequency.one*prefactors.two 
            + mass_frequency.two*prefactors.three
          );
}

//Ansatz for the merger-ringdown amplitude. Equation 19 arXiv:1508.07253:
inline float calculateMergerRingdownAmplitudeAnsatz(
    const useful_powers_s          mass_frequency, 
    const amplitude_coefficients_s coefficients
    ) {
    
    const float ringdown_frequency = coefficients.ringdown_frequency;
    const float damping_time       = coefficients.damping_time;
    
    const float *gamma = coefficients.merger_ringdown;
    
    const float damping_time_gamma3     = damping_time*gamma[2];
    const float fmin_ringdown_frequency = 
        mass_frequency.one - ringdown_frequency;
    
    return expf( -(fmin_ringdown_frequency)*gamma[1] / (damping_time_gamma3) )
        * (damping_time_gamma3*gamma[0]) / (Square(fmin_ringdown_frequency) 
        + Square(damping_time_gamma3));
    }
 
// Ansatz for the intermediate amplitude. Equation 21 arXiv:1508.07253:
inline float calculateIntermediateAmplitudeAnsatz(
    const useful_powers_s          mass_frequency, 
    const amplitude_coefficients_s coefficients
    ) {
    
    const float Mf   = mass_frequency.one;
    const float Mf_2 = Square(Mf);
    
    const float *delta = coefficients.intermediate;
       
    return 
        + delta[0] 
        + Mf*delta[1] 
        + Mf_2*(delta[2] + Mf*delta[3] + Mf_2*delta[4]);
}

// This function computes the IMR amplitude given phenom coefficients. Defined 
// in VIII. Full IMR Waveforms arXiv:1508.07253:
inline float calculateAmplitude(
        const useful_powers_s                  mass_frequency,
        const amplitude_coefficients_s         coefficients,
        const amplitude_inspiral_prefactors_s  prefactors
    ) {
    // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
    // The inspiral, intermediate and merger-ringdown amplitude parts
    // IMRPhenDAmplitude

    const float amplitude_prefactor = 
        prefactors.amp0 * mass_frequency.m_seven_sixths;

    // Split the calculation to just 1 of 3 possible mutually exclusive ranges:
    float amplitude_ansatz = 0.0f;
    
    // Inspiral range:
    if (mass_frequency.one < AMPLITUDE_INTERMEDIATE_START_FREQUENCY) 
    {
        amplitude_ansatz = 
            calculateInspiralAmplitudeAnsatz(
                mass_frequency, prefactors
            );
    }
    // Merger-Ringdown range:
    else if (mass_frequency.one >= coefficients.inspiral_merger_peak_frequency) 
    {
        amplitude_ansatz = 
            calculateMergerRingdownAmplitudeAnsatz(
                mass_frequency, coefficients
            );
    }
    // Intermediate range:
    else  
    {
        amplitude_ansatz = 
            calculateIntermediateAmplitudeAnsatz(
                mass_frequency, coefficients
            );
    }
    
    return amplitude_prefactor*amplitude_ansatz;
}

// Take the AmpInsAnsatz expression and compute the first derivative
// with respect to frequency to get the expression below:
float calculateInspiralAmplitudeAnsatzDerivative(
    useful_powers_s          mass_frequency, 
    amplitude_coefficients_s coefficients
    ) {
    
    const float Mf   = mass_frequency.one;
    const float eta  = coefficients.eta;
    const float chi1 = coefficients.chi1;
    const float chi2 = coefficients.chi2;
    const float rho1 = coefficients.inspiral[0];
    const float rho2 = coefficients.inspiral[1];
    const float rho3 = coefficients.inspiral[2];

    const float chi12     = coefficients.chi12;
    const float chi22     = coefficients.chi22;
    const float eta2      = coefficients.eta2;
    const float eta3      = coefficients.eta3;
    const float Pi        = PI_POWER_ONE;
    const float Pi2       = PI_POWER_TWO;
    const float Seta      = coefficients.Seta;
    const float SetaPlus1 = coefficients.SetaPlus1;

    return 
        (
            (
              - 969.0f 
              + 1804.0f*eta
            )*PI_POWER_TWO_THIRDS
        )/(1008.0f*mass_frequency.third)
        + (
             (chi1*(81.0f*SetaPlus1 - 44.0f*eta) 
            + chi2*(81.0f - 81.0f*Seta - 44.0f*eta))*Pi
        )/48.0f
        + (
            (
                - 27312085.0f 
                - 10287648.0f*chi22 
                - 10287648.0f*chi12*SetaPlus1
                + 10287648.0f*chi22*Seta
                + 24.0f*(
                    - 1975055.0f 
                    + 857304.0f*chi12 
                    - 994896.0f*chi1*chi2 
                    + 857304.0f*chi22
                    )*eta
                + 35371056.0f*eta2
            )*mass_frequency.third*PI_POWER_FOUR_THIRDS
        )/6.096384e6f
        + (
            5.0f*mass_frequency.two_thirds*PI_POWER_FIVE_THIRDS*(chi2*(
                - 285197.0f*(-1.0f + Seta)
                + 4.0f*(-91902.0f + 1579.0f*Seta)*eta 
                - 35632.0f*eta2
            ) 
            + chi1*(
                  285197.0f*SetaPlus1
                - 4.0f*(91902.0f + 1579.0f*Seta)*eta - 35632.0f*eta2
            ) 
            + 42840.0f*(-1.0f + 4.0f*eta)*Pi)
        )/96768.0f
        - (Mf*Pi2*(
                -336.0f*(
                    - 3248849057.0f 
                    + 2943675504.0f*chi12 
                    - 3339284256.0f*chi1*chi2 
                    + 2943675504.0f*chi22
                )*eta2 
                - 324322727232.0f*eta3
                - 7.0f*(
                    - 177520268561.0f 
                    + 107414046432.0f*chi22 
                    + 107414046432.0f*chi12*SetaPlus1 
                    - 107414046432.0f*chi22*Seta
                    + 11087290368.0f*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi
                )
                + 12.0f*eta*(
                    - 545384828789.0f
                    - 176491177632.0f*chi1*chi2 
                    + 202603761360.0f*chi22 
                    + 77616.0f*chi12*(2610335.0f + 995766.0f*Seta)
                    - 77287373856.0f*chi22*Seta 
                    + 5841690624.0f*(chi1 + chi2)*Pi 
                    + 21384760320.0f*Pi2
                )
            )
        )/3.0042980352e10f
        + (7.0f/3.0f)*mass_frequency.four_thirds*rho1 
        + (8.0f/3.0f)*mass_frequency.five_thirds*rho2 
        + 3.0f*mass_frequency.two*rho3;
}


// First frequency derivative of AmpMRDAnsatz:
inline float calculateMergerRingdownAmplitudeAnsatzDerivative(
    const useful_powers_s          mass_frequency,  
    const amplitude_coefficients_s coefficients
    ) {
    
    const float ringdown_frequency = coefficients.ringdown_frequency;
    const float damping_time       = coefficients.damping_time;
    
    const float *gamma = coefficients.merger_ringdown;

    const float damping_time_gamma3      = damping_time * gamma[2];
    const float pow2_damping_time_gamma3 = Square(damping_time_gamma3);
    
    const float fmin_ringdown_frequency  = 
        mass_frequency.one - ringdown_frequency;
    const float expfactor = 
        expf(((fmin_ringdown_frequency)*gamma[1])/(damping_time_gamma3));
    const float pow2pluspow2 = 
        Square(fmin_ringdown_frequency) + pow2_damping_time_gamma3;

    return 
        ((-2.0f*damping_time*(fmin_ringdown_frequency)*gamma[2]*gamma[0]) 
        / pow2pluspow2 - (gamma[1]*gamma[0])) / ( expfactor * (pow2pluspow2));
}


// Formula to predict the final spin. Equation 3.6 arXiv:1508.07250
// s defined around Equation 3.6.
inline float _calcFinalSpin(
    const float eta, 
    const float s
    ) {
    
    const float eta2 = eta*eta;
    const float eta3 = eta2*eta;
    const float s2   = s*s;
    const float s3   = s2*s;

    return eta*(3.4641016151377544f - 4.399247300629289f*eta +
           9.397292189321194f*eta2 - 13.180949901606242f*eta3 +
           s*((1.0f/eta - 0.0850917821418767f - 5.837029316602263f*eta) +
           (0.1014665242971878f - 2.0967746996832157f*eta)*s +
           (-1.3546806617824356f + 4.108962025369336f*eta)*s2 +
           (-0.8676969352555539f + 2.064046835273906f*eta)*s3));
}

inline float calcFinalSpin(
    const system_properties_s system_properties
    ) {
    
    const float eta  = system_properties.symmetric_mass_ratio;
    const float chi1 = system_properties.companion[0].spin.z;
    const float chi2 = system_properties.companion[1].spin.z;
    
    // Convention m1 >= m2
    const float Seta = sqrtf(1.0f - 4.0f*eta);
    const float m1 = 0.5f * (1.0f + Seta);
    const float m2 = 0.5f * (1.0f - Seta);
    const float m1s = m1*m1;
    const float m2s = m2*m2;
    // s defined around Equation 3.6 arXiv:1508.07250:
    const float s = (m1s * chi1 + m2s * chi2);
    
    return _calcFinalSpin(eta, s);
}

useful_powers_s initUsefulPowers(
    const float one
    ) {
        
    const float sixth = sqrtf(one) / cbrtf(one);
    const float m_sixth = 1.0f/sixth;
    
    useful_powers_s powers;
    
    powers.one            = one;
    powers.third          = sixth*sixth;
    powers.two_thirds     = powers.third*powers.third;
    powers.four_thirds    = one*powers.third;
    powers.five_thirds    = powers.four_thirds*powers.third;
    powers.two            = one*one;
    powers.seven_thirds   = powers.third*powers.two;
    powers.eight_thirds   = powers.two_thirds*powers.two;
    
    powers.inv            = 1.0f/one;
    powers.m_seven_sixths = powers.inv*m_sixth;
    powers.m_third        = m_sixth * m_sixth;
    powers.m_two_thirds   = powers.m_third * powers.m_third;
    powers.m_five_thirds  = powers.inv*powers.m_two_thirds;
    
    return powers;
}

float calculateAmplitudePeak(
    const amplitude_coefficients_s amplitude_coefficients
    ) {
    
    const float ringdown_frequency = amplitude_coefficients.ringdown_frequency;
    const float damping_time       = amplitude_coefficients.damping_time;
    
    const float *gamma = amplitude_coefficients.merger_ringdown;

    // NOTE: There's a problem with this expression from the paper becoming 
    // imaginary if gamma2>=1 Fix: if gamma2 >= 1 then set the Square root term 
    // to zero.
    
    float amplitude_peak = 0.0f;
    if (!(gamma[1] > 1.0f))
    {
        amplitude_peak = 
              fabsf(ringdown_frequency 
            + (damping_time*
                (-1 + sqrtf(1 - Square(gamma[1])))*gamma[2])/gamma[1]);
    }
    else
    {
        amplitude_peak = 
            fabsf(ringdown_frequency + (-damping_time*gamma[2])/gamma[1]);
    }
    
    return amplitude_peak;
}

/**************************** Amplitude functions *****************************/
 

// Amplitude scaling factor defined by eq. 17 in 1508.07253:
inline float calculateAmplitude0(
    const float symmetric_mass_ratio
    ) {
    return sqrtf(2.0f/3.0f*symmetric_mass_ratio)*PI_M_SIXTH;
}

amplitude_inspiral_prefactors_s initAmplitudeInspiralPrefactors(
    const amplitude_coefficients_s coefficients
    ) {
    
    amplitude_inspiral_prefactors_s prefactors;

    const float eta = coefficients.eta;

    prefactors.amp0 = calculateAmplitude0(eta);

    const float chi1 = coefficients.chi1;
    const float chi2 = coefficients.chi2;

    const float chi12 = coefficients.chi12;
    const float chi22 = coefficients.chi22;
    const float eta2  = coefficients.eta2;
    const float eta3  = coefficients.eta3;

    const float Pi        = (float)M_PI;
    const float Pi2       = PI_POWER_TWO;
    const float Seta      = coefficients.Seta;
    const float SetaPlus1 = coefficients.SetaPlus1;

    prefactors.two_thirds = 
        ((-969.0f + 1804.0f*eta)*PI_POWER_TWO_THIRDS)/672.0f;
    
    prefactors.one = 
        (
             (chi1*(81.0f*SetaPlus1 - 44.0f*eta) 
            + chi2*(81.0f - 81.0f*Seta - 44.0f*eta))*Pi
        ) 
        / 48.0f;
        
    prefactors.four_thirds = 
        (
            (
                - 27312085.0f 
                - 10287648.0f*chi22 
                - 10287648.0f*chi12*SetaPlus1 
                + 10287648.0f*chi22*Seta
                + 24.0f*(
                    - 1975055.0f 
                    + 857304.0f*chi12 
                    - 994896.0f*chi1*chi2 
                    + 857304.0f*chi22
                    )*eta
                + 35371056.0f*eta2
            ) 
            * PI_POWER_FOUR_THIRDS
        ) 
        / 8.128512e6f;
        
    prefactors.five_thirds = 
        (
            PI_POWER_FIVE_THIRDS*(
                chi2*(
                    - 285197.0f*(-1.0f + Seta) 
                    + 4.0f*(-91902.0f + 1579.0f*Seta)*eta 
                    - 35632.0f*eta2
                    )
                + chi1*(
                      285197.0f*SetaPlus1 
                    - 4.0f*(91902.0f + 1579.0f*Seta)*eta 
                    - 35632.0f*eta2
                    )
                + 42840.0f*(-1.0f + 4.0f*eta)*Pi
            )
        ) / 32256.0f;
        
    prefactors.two = 
        -1.0f*(
            Pi2*(
                - 336.0f*(
                    - 3248849057.0f 
                    + 2943675504.0f*chi12 
                    - 3339284256.0f*chi1*chi2 
                    + 2943675504.0f*chi22
                    )*eta2
                - 324322727232.0f*eta3
                - 7.0f*(
                    - 177520268561.0f 
                    + 107414046432.0f*chi22 
                    + 107414046432.0f*chi12*SetaPlus1
                    - 107414046432.0f*chi22*Seta 
                    + 11087290368.0f*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi
                    )
                + 12.0f*eta*(
                    - 545384828789.0f 
                    - 176491177632.0f*chi1*chi2 
                    + 202603761360.0f*chi22
                    + 77616.0f*chi12*(2610335.0f + 995766.0f*Seta) 
                    - 77287373856.0f*chi22*Seta
                    + 5841690624.0f*(chi1 + chi2)*Pi 
                    + 21384760320.0f*Pi2
                 )
            )
        ) / 6.0085960704e10f;
        
    prefactors.seven_thirds = coefficients.inspiral[0];
    prefactors.eight_thirds = coefficients.inspiral[1];
    prefactors.three        = coefficients.inspiral[2];

    return prefactors;
}

amplitude_coefficients_s computeDeltasFromCollocation(
    amplitude_coefficients_s amplitude_coefficients
    ) {

    const int32_t num_collocation_points = 3;
    // Three evenly spaced collocation points in the interval [f1,f3]:
    
    const float collocation_mid_point = 
        0.5f*
        (amplitude_coefficients.inspiral_merger_peak_frequency 
         - AMPLITUDE_INTERMEDIATE_START_FREQUENCY);
    
    const float collocation_points[] = 
    {
        AMPLITUDE_INTERMEDIATE_START_FREQUENCY,
        collocation_mid_point,
        amplitude_coefficients.inspiral_merger_peak_frequency
    };
    
    useful_powers_s collocation_point_powers[num_collocation_points];
    for (int32_t index = 0; index < num_collocation_points; index++)
    {
        collocation_point_powers[index] = 
            initUsefulPowers(collocation_points[index]);    
    }

    amplitude_inspiral_prefactors_s prefactors = 
        initAmplitudeInspiralPrefactors(amplitude_coefficients);

    // v1 is inspiral model evaluated at f1
    // d1 is derivative of inspiral model evaluated at f1.
    const float v1 = 
        calculateInspiralAmplitudeAnsatz(
            collocation_point_powers[0], 
            prefactors
        );
    
    const float d1 = 
        calculateInspiralAmplitudeAnsatzDerivative(
            collocation_point_powers[0], 
            amplitude_coefficients
        );
    
    // v3 is merger-ringdown model evaluated at f3
    // d2 is derivative of merger-ringdown model evaluated at f3.
    const float v3 = 
        calculateMergerRingdownAmplitudeAnsatz(
            collocation_point_powers[2], 
            amplitude_coefficients
        );
    const float d2 = 
        calculateMergerRingdownAmplitudeAnsatzDerivative(
            collocation_point_powers[2], 
            amplitude_coefficients
        );

    // v2 is the value of the amplitude evaluated at f2
    // they come from the fit of the collocation points in the intermediate 
    // region.
    const float v2 = 
       calcCoefficient(
           AMPLITUDE_INTERMEDIATE_COLLOCATION_FIT_COEFFICIENT_TERMS, 
           amplitude_coefficients.eta, 
           amplitude_coefficients.eta2, 
           (-1.0f + amplitude_coefficients.chi)
       );  
    
    // Now compute the delta_i's from the collocation coefficients
    // Precompute common quantities here and pass along to delta functions:
    const int32_t num_powers = 6;
    float f[num_collocation_points][num_powers];
    for (int32_t index = 0; index < num_collocation_points; index++)
    {
        f[index][1] = collocation_points[index];
        f[index][2] = f[index][1] * f[index][1];
        f[index][3] = f[index][2] * f[index][1];
        f[index][4] = f[index][3] * f[index][1];
        f[index][5] = f[index][4] * f[index][1];
    }
    
    // Use inner term macros defined in phenomd_data.h to abstract maths:
    const float inner_term_results[] =
    {
        NUM_AMPLITUDE_INTERMEDIATE_TERMS_0,
        NUM_AMPLITUDE_INTERMEDIATE_TERMS_1,
        NUM_AMPLITUDE_INTERMEDIATE_TERMS_2,
        NUM_AMPLITUDE_INTERMEDIATE_TERMS_3,
        NUM_AMPLITUDE_INTERMEDIATE_TERMS_4
    };
    
    const float divisor =
          Square(f[0][1] - f[1][1])
        * Cube  (f[0][1] - f[2][1])
        * Square(f[2][1] - f[1][1]);
    
    for (int32_t index = 0; index < NUM_AMPLITUDE_INTERMEDIATE_TERMS; index++)
    {
        amplitude_coefficients.intermediate[index] 
            = -1.0f*inner_term_results[index] / divisor;
    }     
    
    return amplitude_coefficients;
}

amplitude_coefficients_s initAmplitudeCoefficients(
    const system_properties_s system_properties,
    const float               finspin
    ) {
    
    // Unpack properties:
    const float eta     = system_properties.symmetric_mass_ratio;
    const float chi1    = system_properties.companion[0].spin.z;
    const float chi2    = system_properties.companion[1].spin.z;
    
    amplitude_coefficients_s coefficients;

    coefficients.eta = eta;
    coefficients.etaInv = 1.0f/eta;
    coefficients.chi1 = chi1;
    coefficients.chi2 = chi2;
    coefficients.chi12 = chi1*chi1;
    coefficients.chi22 = chi2*chi2;
    const float eta2 = eta*eta;
    coefficients.eta2 = eta2;
    coefficients.eta3 = eta*eta2;
    const float Seta = sqrtf(1.0f - 4.0f*eta);
    coefficients.Seta = Seta;
    coefficients.SetaPlus1 = 1.0f + Seta;

    coefficients.q = 0.5f * (1.0f + Seta - 2.0f*eta) * coefficients.etaInv;
    coefficients.chi = chiPN(Seta, eta, chi1, chi2);
    const float xi = -1.0f + coefficients.chi;

    coefficients.ringdown_frequency = 
        calculateRingdownFrequency(eta, chi1, chi2, finspin, QNMData_fring);
    coefficients.damping_time = 
        calculateRingdownFrequency(eta, chi1, chi2, finspin, QNMData_fdamp);
    
    // Compute gamma_i's, rho_i's first then delta_i's:
    calcCoefficients(
        AMPLITUDE_MERGER_RINGDOWN_COEFFICIENT_TERMS,
        NUM_AMPLITUDE_MERGER_RINGDOWN_COEFFICIENTS,
        eta, eta2, xi,
        coefficients.merger_ringdown
    );
    calcCoefficients(
        AMPLITUDE_INSPIRAL_COEFFICIENT_TERMS,
        NUM_AMPLITUDE_INSPIRAL_COEFFICIENTS,
        eta, eta2, xi,
        coefficients.inspiral
    );

    coefficients.inspiral_merger_peak_frequency 
        = calculateAmplitudePeak(coefficients);

    // Compute delta_is:
    coefficients = 
        computeDeltasFromCollocation(coefficients);
        
    return coefficients;
}


// Ansatz for the inspiral phase. We call the LAL TF2 coefficients here. The 
// exact values of the coefficients used are given as comments in the top of 
// this file. Defined by Equation 27 and 28 arXiv:1508.07253
inline float calculateInspiralPhaseAnsatz(
    const useful_powers_s             mass_frequency, 
    const phase_inspiral_prefactors_s prefactors, 
    const phase_coefficients_s        coefficients
    ) {
    
    // Assemble PN phasing series:
    const float v    = mass_frequency.third * CUBE_ROOT_PI;
    const float logv = logf(v);
    
    const float phase_ansatz =  
         prefactors.initial_phasing
       + prefactors.two_thirds        * mass_frequency.two_thirds
       + prefactors.third             * mass_frequency.third
       + prefactors.third_with_logv   * logv * mass_frequency.third
       + prefactors.logv              * logv
       + prefactors.minus_third       * mass_frequency.m_third
       + prefactors.minus_two_thirds  * mass_frequency.m_two_thirds
       + prefactors.minus_one         * mass_frequency.inv
       + prefactors.minus_four_thirds / mass_frequency.four_thirds
       + prefactors.minus_five_thirds * mass_frequency.m_five_thirds
       // Now add higher order terms that were calibrated for PhenomD:
       + coefficients.etaInv*
            (    prefactors.one         * mass_frequency.one 
               + prefactors.four_thirds * mass_frequency.four_thirds
               + prefactors.five_thirds * mass_frequency.five_thirds
               + prefactors.two         * mass_frequency.two
            );
        
    return phase_ansatz;
}

inline float calculateIntermediatePhaseAnsatzDerivative(
    const useful_powers_s      mass_frequency, 
    const phase_coefficients_s coefficients
    ) {
    
    const float *beta = coefficients.intermediate;
    
   return 
       (
             beta[0] 
           + beta[2]/Quart(mass_frequency.one) 
           + beta[1]/mass_frequency.one
        )*coefficients.etaInv;
}

// Ansatz for the intermediate phase defined by Equation 16 arXiv:1508.07253
inline float calculateIntermediatePhaseAnsatz(
    const useful_powers_s      mass_frequency, 
    const phase_coefficients_s coefficients
    ) {
    // 1./eta in paper omitted and put in when need in the functions:
    // ComputeIMRPhenDPhaseConnectionCoefficients
    // IMRPhenDPhase
    
    const float *beta = coefficients.intermediate;
    
    return 
          beta[0]*mass_frequency.one 
        - beta[2]/(3.0f*Cube(mass_frequency.one)) 
        + beta[1]*logf(mass_frequency.one);
}

//Ansatz for the merger-ringdown phase Equation 14 arXiv:1508.07253
//Rholm was added when IMRPhenomHM (high mode) was added.
//Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
//to 1. and PhenomD is recovered.
//Taulm = damping_timelm/damping_time22. Ratio of ringdown damping times.
//Again, when Taulm = 1.0 then PhenomD is recovered.
inline float calculateMergerRingdownPhaseAnsatz(
    const useful_powers_s      mass_frequency, 
    const phase_coefficients_s coefficients, 
    const float                Rholm, 
    const float                Taulm
    ) {
    
    const float sqrootf  = sqrtf(mass_frequency.one);
    const float fpow1_5  = mass_frequency.one * sqrootf;
    const float fpow0_75 = sqrtf(fpow1_5);
    
    const float *alpha = coefficients.merger_ringdown;
    
    return 
        - (alpha[1]/mass_frequency.one)
        + (4.0f/3.0f) * (alpha[2] * fpow0_75)
        + alpha[0] * mass_frequency.one
        + alpha[3] * Rholm * atanf(
                  (mass_frequency.one - alpha[4] 
                  * coefficients.ringdown_frequency) 
                / (Rholm * coefficients.damping_time * Taulm)
            );
}


// This function computes the IMR phase given phenom coefficients.
// Defined in VIII. Full IMR Waveforms arXiv:1508.07253
// Rholm was added when IMRPhenomHM (high mode) was added.
// Rholm = ringdown_frequency22/ringdown_frequencylm. For PhenomD (only 
// (l,m)=(2,2)) this is just equal to 1. and PhenomD is recovered.
// Taulm = damping_timelm/damping_time22. Ratio of ringdown damping times.
// Again, when Taulm = 1.0 then PhenomD is recovered.
float calculatePhase(
    const useful_powers_s             mass_frequency, 
    const phase_coefficients_s        coefficients, 
    const phase_inspiral_prefactors_s prefactors, 
    const float                       Rholm, 
    const float                       Taulm
    ) {
    
    // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
    // The inspiral, intermendiate and merger-ringdown phase parts

    // Split the calculation to just 1 of 3 possible mutually exclusive ranges
    float phase = 0.0f;
    // Inspiral range
    if (mass_frequency.one < coefficients.intermediate_start_frequency)        
    {
        phase = 
            calculateInspiralPhaseAnsatz(
               mass_frequency, 
               prefactors, 
               coefficients            
            );
    }
    // Merger-Ringdown range:
    else if (mass_frequency.one >= coefficients.merger_ringdown_start_frequency) 
    {
        phase = 
            coefficients.etaInv*
            calculateMergerRingdownPhaseAnsatz(
                mass_frequency, 
                coefficients, 
                Rholm, 
                Taulm
            ) 
            + coefficients.C1MRD 
            + coefficients.C2MRD*mass_frequency.one;
    }
    // Intermediate range:
    else
    {
        phase = 
            coefficients.etaInv*
            calculateIntermediatePhaseAnsatz(
                mass_frequency, 
                coefficients
            ) 
            + coefficients.C1Int 
            + coefficients.C2Int*mass_frequency.one;
    }
    
    return phase;
}

phase_coefficients_s initPhaseCoefficients(
    const system_properties_s system_properties,
    const float               finspin
    ) {
    
    // Unpack properties:
    const float eta  = system_properties.symmetric_mass_ratio;
    const float chi1 = system_properties.companion[0].spin.z;
    const float chi2 = system_properties.companion[1].spin.z;
    
    phase_coefficients_s coefficients;

    // Convention m1 >= m2:
    coefficients.eta = eta;
    coefficients.etaInv = 1.0f/eta;
    coefficients.chi1 = chi1;
    coefficients.chi2 = chi2;
    const float eta2 = eta*eta;
    coefficients.eta2 = eta2;
    coefficients.Seta = sqrtf(1.0f - 4.0f*eta);

    coefficients.q = 
        0.5f*(1.0f + coefficients.Seta - 2.0f*eta)*coefficients.etaInv;
    coefficients.chi = chiPN(coefficients.Seta, eta, chi1, chi2);
    const float xi = -1.0f + coefficients.chi;
    
    calcCoefficients(
        PHASE_INSPIRAL_COEFFICIENT_TERMS,
        NUM_PHASE_INSPIRAL_COEFFICIENTS,
        eta, eta2, xi,
        coefficients.inspiral
    );
    
    calcCoefficients(
        PHASE_INTERMEDIATE_COEFFICIENT_TERMS,
        NUM_PHASE_INTERMEDIATE_COEFFICIENTS,
        eta, eta2, xi,
        coefficients.intermediate
    );
    
    calcCoefficients(
        PHASE_MERGER_RINGDOWN_COEFFICIENT_TERMS,
        NUM_PHASE_MERGER_RINGDOWN_COEFFICIENTS,
        eta, eta2, xi,
        coefficients.merger_ringdown
    );

    coefficients.ringdown_frequency = 
       calculateRingdownFrequency(eta, chi1, chi2, finspin, QNMData_fring);
    coefficients.damping_time = 
       calculateRingdownFrequency(eta, chi1, chi2, finspin, QNMData_fdamp);
       
    return coefficients;
}

phase_inspiral_prefactors_s initPhaseInspiralPrefactors(
    const phase_coefficients_s coefficients,
    const pn_phasing_series_s  series
) {  

    phase_inspiral_prefactors_s prefactors;
    
    const float *sigma = coefficients.inspiral; 
     
    // PN phasing series
    prefactors.initial_phasing   = series.v[5]     - (float)M_PI_4;
    prefactors.two_thirds        = series.v[7]     * PI_POWER_TWO_THIRDS;
    prefactors.third             = series.v[6]     * PI_POWER_ONE_THIRD;
    prefactors.third_with_logv   = series.vlogv[6] * PI_POWER_ONE_THIRD;
    prefactors.logv              = series.vlogv[5];
    prefactors.minus_third       = series.v[4]     * PI_POWER_MINUS_ONE_THIRD;
    prefactors.minus_two_thirds  = series.v[3]     * PI_POWER_MINUS_TWO_THIRDS;
    prefactors.minus_one         = series.v[2]     * PI_POWER_MINUS_ONE;
    prefactors.minus_four_thirds = series.v[1]     / PI_POWER_FOUR_THIRDS;
    prefactors.minus_five_thirds = series.v[0]     * PI_POWER_MINUS_FIVE_THIRDS;
  
    // higher order terms that were calibrated for PhenomD
    prefactors.one         = sigma[0];
    prefactors.four_thirds = sigma[1] * 0.75f;
    prefactors.five_thirds = sigma[2] * 0.6f;
    prefactors.two         = sigma[3] * 0.5f;

    return prefactors;
}

// Temporary instance of DPhiIntAnsatz used when computing coefficients to make 
// the phase C(1) continuous between regions:
inline float calculateTempIntermediatelPhaseAnsatzDerivitive(
    const phase_coefficients_s coefficients
    ) {
    
    const float  start_frequency = coefficients.merger_ringdown_start_frequency;
    const float  etaInv          = coefficients.etaInv;
    const float  C2Int           = coefficients.C2Int;
    const float *beta            = coefficients.intermediate;

    return 
        + C2Int 
        + (
            + beta[0] 
            + beta[2]/Quart(start_frequency) 
            + beta[1]/start_frequency
        )*etaInv;
}

float calculateInspiralPhaseAnsatzDerivitive(
    const phase_coefficients_s coefficients, 
    const pn_phasing_series_s  pn_series
    ) {
    // First frequency derivative of PhiInsAnsatzInt
    
    const float start_freqeuncy = coefficients.intermediate_start_frequency;
    const float *sigma          = coefficients.inspiral; 

    // Assemble PN phasing series:
    const float init_v = cbrtf((float)M_PI*start_freqeuncy);
    
    const int32_t PN_PHASING_SERIES_LENGTH = 9;
    float v[PN_PHASING_SERIES_LENGTH];
    v[0] = 1.0f;
    for (int32_t index = 1; index < PN_PHASING_SERIES_LENGTH + 1; index++)
    {
        v[index] = init_v * v[index-1];
    }
    
    const float logv = logf(init_v);

    // Apply the correct prefactors to LAL phase coefficients to get the
    // phase derivative dphi / dMf = dphi/dv * dv/dMf:
    float phase_ansatz_derivitive = 
        + 2.0f*pn_series.v[7]*v[7]
        + (pn_series.v[6] + pn_series.vlogv[6] * (1.0f + logv)) * v[6]
        + pn_series.vlogv[5] *v[5]
        - 1.0f*pn_series.v[4]*v[4]
        - 2.0f*pn_series.v[3]*v[3]
        - 3.0f*pn_series.v[2]*v[2]
        - 4.0f*pn_series.v[1]*v[1]
        - 5.0f*pn_series.v[0]*v[0];

    phase_ansatz_derivitive /= v[8] * 3.0f;
    phase_ansatz_derivitive *= (float)M_PI;

    // Now add higher order terms that were calibrated for PhenomD
    phase_ansatz_derivitive += 
        (
               sigma[0]*v[0]
             + sigma[1]*v[1]*PI_POWER_MINUS_ONE_THIRD
             + sigma[2]*v[2]*PI_POWER_MINUS_TWO_THIRDS
             + sigma[3]*v[3]*PI_POWER_MINUS_ONE
         )*coefficients.etaInv;

   return phase_ansatz_derivitive;
}


// First frequency derivative of PhiMRDAnsatzInt
// Rholm was added when IMRPhenomHM (high mode) was added.
// Rholm = ringdown_frequency22/ringdown_frequencylm. For PhenomD (only (l,m)=(2,2)) this is just equal
// to 1. and PhenomD is recovered.
// Taulm = damping_timelm/damping_time22. Ratio of ringdown damping times.
// Again, when Taulm = 1.0 then PhenomD is recovered.
float calculateMergerRingdownPhaseAnsatzDerivitive(
    const float                start_frequency, 
    const phase_coefficients_s coefficients, 
    const float                Rholm, 
    const float                Taulm
    ) {

    const float *alpha             = coefficients.merger_ringdown;
    const float ringdown_frequency = coefficients.ringdown_frequency;
    const float damping_time       = coefficients.damping_time;
    
    return (
        + alpha[0] 
        + alpha[1]/Square(start_frequency) 
        + alpha[2]/powf(start_frequency, 0.25f)
        + alpha[3]/(
            damping_time*Taulm*(
                + 1.0f 
                + Square(
                    + start_frequency 
                    - alpha[4]* 
                    ringdown_frequency
                )/Square(
                    damping_time*Taulm*Rholm
                )
            )
        ) 
    )*coefficients.etaInv;
}

phase_coefficients_s computePhenomDPhaseConnectionCoefficients(
    phase_coefficients_s         coefficients, 
    pn_phasing_series_s          pn_series, 
    phase_inspiral_prefactors_s  prefactors, 
    const float Rholm, 
    const float Taulm
) {
    
    const float etaInv = coefficients.etaInv;

    // Transition frequencies
    // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
    coefficients.intermediate_start_frequency = 
        PHASE_INTERMEDIATE_START_FREQUENCY;
    coefficients.merger_ringdown_start_frequency = 
        0.5f*coefficients.ringdown_frequency;

    // Compute C1Int and C2Int coeffs
    // Equations to solve for to get C(1) continuous join
    // PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
    // Joining at fInsJoin
    // PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
    // PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
    
    const useful_powers_s intermediate_start_frequency_powers = 
        initUsefulPowers(coefficients.intermediate_start_frequency);
    const useful_powers_s merger_ringdown_start_frequency_powers = 
        initUsefulPowers(coefficients.merger_ringdown_start_frequency);
        
    const float inspiral_ansatz_derivitive = 
        calculateInspiralPhaseAnsatzDerivitive(
            coefficients, 
            pn_series
        );
        
    const float intermediate_ansatz_derivitive = 
        calculateIntermediatePhaseAnsatzDerivative(
            intermediate_start_frequency_powers, 
            coefficients
        );
    
    coefficients.C2Int = 
        inspiral_ansatz_derivitive - intermediate_ansatz_derivitive;
    
    const float inspiral_ansatz = 
        calculateInspiralPhaseAnsatz(
            intermediate_start_frequency_powers, 
            prefactors, 
            coefficients
        );
        
    const float intermediate_ansatz = 
        calculateIntermediatePhaseAnsatz(
            intermediate_start_frequency_powers, 
            coefficients
        );
    
    coefficients.C1Int = 
        + inspiral_ansatz
        - etaInv*intermediate_ansatz
        - coefficients.C2Int*coefficients.intermediate_start_frequency;

    // Compute C1MRD and C2MRD coeffs
    // Equations to solve for to get C(1) continuous join
    // PhiInsInt (f)  =   PhiMRD (f) + C1MRD + C2MRD f
    // Joining at merger_ringdown_start_frequency
    // Where \[Phi]InsInt(f) is the \[Phi]Ins+\[Phi]Int joined function
    // PhiInsInt (merger_ringdown_start_frequency) 
    //   = PhiMRD (merger_ringdown_start_frequency) 
    //   + C1MRD + C2MRD merger_ringdown_start_frequency
    // PhiInsInt'(merger_ringdown_start_frequency)  
    //  = PhiMRD'(merger_ringdown_start_frequency) + C2MRD
    // temporary Intermediate Phase function to Join up the Merger-Ringdown
    const float temp_intermediate_ansatz = 
        etaInv*calculateIntermediatePhaseAnsatz(
            merger_ringdown_start_frequency_powers, 
            coefficients
        ) 
        + coefficients.C1Int 
        + coefficients.C2Int*coefficients.merger_ringdown_start_frequency;
        
    const float temp_intermediate_ansatz_derivitive =
        calculateTempIntermediatelPhaseAnsatzDerivitive(
            coefficients
        );
        
    const float merger_ringdown_ansatz_derivitive = 
        calculateMergerRingdownPhaseAnsatzDerivitive(
            coefficients.merger_ringdown_start_frequency,
            coefficients, 
            Rholm, 
            Taulm
        );

    coefficients.C2MRD = 
        + temp_intermediate_ansatz_derivitive 
        - merger_ringdown_ansatz_derivitive;
    coefficients.C1MRD = 
        temp_intermediate_ansatz 
        - etaInv*calculateMergerRingdownPhaseAnsatz(
            merger_ringdown_start_frequency_powers, 
            coefficients, 
            Rholm, 
            Taulm) 
        - coefficients.C2MRD*coefficients.merger_ringdown_start_frequency;
        
    return coefficients;
}