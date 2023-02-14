#ifndef _PHENOM_H
#define _PHENOM_H

#include <omp.h>

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#define LAL_PI 3.141592653589793238462643383279502884
#define LAL_MTSUN_SI 4.925491025543575903411922162094833998e-6
#define LAL_MSUN_SI 1.988409902147041637325262574352366540e30
#define LAL_MRSUN_SI 1.476625061404649406193430731479084713e3

// Fraction of waveform duration to add as extra time for tapering:
#define EXTRA_TIME_FRACTION 0.1
// More extra time measured in cycles at the starting frequency:
#define EXTRA_CYCLES 3.0 

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
     const int32_t             num_strain_axis_samples,
     const frequencyUnit_t     starting_frequency,
     const frequencyUnit_t     ending_frequency,
     const frequencyUnit_t     frequency_interval,      
     const angularUnit_t       reference_phase,      // phase at reference frequency
     const frequencyUnit_t     reference_frequency,  // reference frequency 
     const system_properties_s system_properties,
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
    
    const float m1 = system_properties.companion[0].mass.msun;
    const float m2 = system_properties.companion[1].mass.msun;

    const float chi1 = system_properties.companion[0].spin.z;
    const float chi2 = system_properties.companion[1].spin.z;
    
    const float distance = system_properties.distance.meters;
    const float inclination = system_properties.inclination.radians;
            
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
        offset = (size_t) (starting_frequency.hertz / frequency_interval.hertz);
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
            final_spin
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
    
    // Compute coefficients to make phase C^1 continuous (phase and first derivative)
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

    // incorporating reference frequency
    const float reference_mass_frequency = total_mass.seconds * reference_frequency.hertz;
    useful_powers_s powers_of_fRef = initUsefulPowers(reference_mass_frequency);
    
    const float phifRef = 
        calculatePhase(
            powers_of_fRef, 
            phase_coefficients, 
            phase_prefactors,
            Rholm, 
            Taulm
        );

    // factor of 2 b/c reference_phase is orbital phase
    const float phi_precalc = 2.0f*reference_phase.radians + phifRef;
    
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
    const float phi0,                  /**< Orbital phase at fRef (rad) */
    const float fRef_in,               /**< reference frequency (Hz) */
    const float deltaF,                /**< Sampling frequency (Hz) */
    const float m1_SI,                 /**< Mass of companion 1 (kg) */
    const float m2_SI,                 /**< Mass of companion 2 (kg) */
    const float chi1,                  /**< Aligned-spin parameter of companion 1 */
    const float chi2,                  /**< Aligned-spin parameter of companion 2 */
    const float f_min,                 /**< Starting GW frequency (Hz) */
    const float f_max,                 /**< End frequency; 0 defaults to Mf = \ref f_CUT */
    const float distance,               /**< Distance of source (m) */
    const float inclination,
    LALDict *extraParams, /**< linked list containing the extra testing GR parameters */
    NRTidal_version_type  NRTidal_version /**< Version of NRTides; can be one of NRTidal versions or NoNRT_V for the BBH baseline */
    ) {
    
    /* external: SI; internal: solar masses */
    const float m1 = m1_SI / LAL_MSUN_SI;
    const float m2 = m2_SI / LAL_MSUN_SI;
    
    /* check inputs for sanity */
    //if (fRef_in < 0) //XLAL_ERROR(XLAL_EDOM, "fRef_in must be positive (or 0 for 'ignore')\n");
    //if (deltaF <= 0) //XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
    //if (m1 <= 0) //XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
    //if (m2 <= 0) //XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");
    // if (f_min <= 0) //XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
    // if (f_max < 0) //XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");
    // if (distance <= 0) //XLAL_ERROR(XLAL_EDOM, "distance must be positive\n");

    const float q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

    if (q > 1000)
     XLAL_PRINT_WARNING("Warning: The model is not supported for high mass ratio, see MAX_ALLOWED_MASS_RATIO\n");

   // if (chi1 > 1.0 || chi1 < -1.0 || chi2 > 1.0 || chi2 < -1.0)
     //XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    // if no reference frequency given, set it to the starting GW frequency
    float fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    const float M_sec = (m1+m2) * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
    const float fCut = f_CUT/M_sec; // convert Mf -> Hz
    // Somewhat arbitrary end point for the waveform.
    // Chosen so that the end of the waveform is well after the ringdown.
    //if (fCut <= f_min)
     //XLAL_ERROR(XLAL_EDOM, "(fCut = %g Hz) <= f_min = %g\n", fCut, f_min);

     /* default f_max to Cut */
    float f_max_prime = f_max;
    f_max_prime = f_max ? f_max : fCut;
    f_max_prime = (f_max_prime > fCut) ? fCut : f_max_prime;
    //if (f_max_prime <= f_min)
     //XLAL_ERROR(XLAL_EDOM, "f_max <= f_min\n");
    
    int32_t num_strain_axis_samples = 0; timeUnit_t gps_time = initTimeSeconds(0.0f);

    // Set up output array with size closest power of 2:
    num_strain_axis_samples = calcNextPow2(f_max / deltaF) + 1;

    // Coalesce at gps_time = 0:
    // Shift by overall length in time:  
    gps_time = addTimes(2, gps_time, initTimeSeconds(-1. / deltaF));
    
    if (f_max_prime < f_max) 
    {
         // The user has requested a higher f_max than Mf=fCut.
         // Resize the frequency series to fill with zeros beyond the cutoff frequency.
        num_strain_axis_samples = calcNextPow2(f_max / deltaF) + 1; // we actually want to have the length be a power of 2 + 1
    }
    
    spin_t spin_a = {
        .x = 0,
        .y = 0,
        .z = chi1
    };
    
    spin_t spin_b = {
        .x = 0,
        .y = 0,
        .z = chi2
    };
    
    companion_s companion_a = 
    {
        .mass              = initMassKilograms(m1_SI),
        .spin              = spin_a,
        .quadrapole_moment = 0.0,
        .lambda            = 0.0
    };
    companion_s companion_b = 
    {
        .mass              = initMassKilograms(m2_SI),
        .spin              = spin_b,
        .quadrapole_moment = 0.0,
        .lambda            = 0.0
    };
    
    system_properties_s system_properties =
        initBinarySystem(
            companion_a,
            companion_b,
            initLengthMeters(distance),
            0.0,
            initAngleRadians(inclination),
            initAngleRadians(0.0),
            0.0,
            0.0, 
            0.0
        );
    
    complex_waveform_axes_s waveform_axes = 
        _cuPhenomDGenerateFD(
            num_strain_axis_samples,
            initFrequencyHertz(f_min),
            initFrequencyHertz(f_max_prime),
            initFrequencyHertz(deltaF), 
            initAngleRadians(phi0),
            initFrequencyHertz(fRef),
            system_properties,
            NRTidal_version
        );
    
    return waveform_axes;
}

int _convertWaveformAxesToLAL_temp(
    complex_waveform_axes_s waveform_axes,
    float deltaF,
    float f_min,
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde     /**< FD cross polarization */
    ) {
    
    /* Produce both polarizations */
    const LALUnit lalStrainUnit = { 0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
    
    if (deltaF > 0) // Freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0:
    { 
        *hptilde = 
            XLALCreateCOMPLEX16FrequencySeries(
                "htilde: FD waveform", 
                (LIGOTimeGPS){ 0, 0 }, 
                0.0, 
                deltaF, 
                &lalStrainUnit, 
                waveform_axes.strain.num_samples
            );

        *hctilde = 
            XLALCreateCOMPLEX16FrequencySeries(
                "ctilde: FD waveform", 
                (LIGOTimeGPS){ 0, 0 }, 
                0.0, 
                deltaF, 
                &lalStrainUnit, 
                waveform_axes.strain.num_samples
            );
    }
    else // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
    {
        *hptilde = 
            XLALCreateCOMPLEX16FrequencySeries(
                "htilde: FD waveform", 
                (LIGOTimeGPS){ 0, 0 },
                f_min, 
                deltaF, 
                &lalStrainUnit, 
                waveform_axes.strain.num_samples
            );
        *hctilde = 
            XLALCreateCOMPLEX16FrequencySeries(
                "ctilde: FD waveform", 
                (LIGOTimeGPS){ 0, 0 }, 
                f_min, 
                deltaF, 
                &lalStrainUnit, 
                waveform_axes.strain.num_samples
            );
    }

    complex_strain_element_c *strain_fd = NULL;
    cudaToHost(
        waveform_axes.strain.values, 
        sizeof(complex_strain_element_c),
        waveform_axes.strain.num_samples,
        &strain_fd
    );


    for(int32_t j = 0; j < (*hptilde)->data->length; j++) 
    {    
     (*hctilde)->data->data[j] = (complex double)strain_fd[j].cross;
     (*hptilde)->data->data[j] = (complex double)strain_fd[j].plus;
    }
}

int cuInspiralChooseFDWaveform(
     COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
     COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
     const double m1,                         /**< mass of companion 1 (kg) */
     const double m2,                         /**< mass of companion 2 (kg) */
     const double S1x,                        /**< x-component of the dimensionless spin of object 1 */
     const double S1y,                        /**< y-component of the dimensionless spin of object 1 */
     const double S1z,                        /**< z-component of the dimensionless spin of object 1 */
     const double S2x,                        /**< x-component of the dimensionless spin of object 2 */
     const double S2y,                        /**< y-component of the dimensionless spin of object 2 */
     const double S2z,                        /**< z-component of the dimensionless spin of object 2 */
     const double distance,                   /**< distance of source (m) */
     const double inclination,                /**< inclination of source (rad) */
     const double phiRef,                     /**< reference orbital phase (rad) */
     const double longAscNodes,               /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
     const double eccentricity,               /**< eccentricity at reference epoch */
     const double meanPerAno,                 /**< mean anomaly of periastron */
     // frequency sampling parameters, no default value
     const double deltaF,                     /**< sampling interval (Hz) */
     const double f_min,                      /**< starting GW frequency (Hz) */
     const double f_max,                      /**< ending GW frequency (Hz) */
     double f_ref,                            /**< Reference frequency (Hz) */
     LALDict *LALparams,                     /**< LAL dictionary containing accessory parameters */
     const Approximant approximant           /**< post-Newtonian approximant to use for waveform production */
     ) {
     int ret;
     unsigned int j;
     double pfac, cfac;
  
     /* The non-precessing waveforms return h(f) for optimal orientation
     * (i=0, Fp=1, Fc=0; Lhat pointed toward the observer)
     * To get generic polarizations we multiply by inclination dependence
     * and note hc(f) \propto -I * hp(f)
     * Non-precessing waveforms multiply hp by pfac, hc by -I*cfac
     */
     cfac = cos(inclination);
     pfac = 0.5 * (1. + cfac*cfac);
    
    complex float *strain_fd_g = NULL;
    
    complex_waveform_axes_s waveform_axes;
  
     switch (approximant)
     {
         case IMRPhenomD:
             // Call the waveform driver routine:
            waveform_axes = 
                 cuPhenomDGenerateFD(
                     phiRef, 
                     f_ref, 
                     deltaF, 
                     m1, 
                     m2,
                     S1z, 
                     S2z, 
                     f_min, 
                     f_max, 
                     distance, 
                     inclination,
                     LALparams, 
                     4
                 );
             
             const LALUnit lalStrainUnit = { 0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
             if (deltaF > 0) // Freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0:
                { 
                    *hptilde = 
                        XLALCreateCOMPLEX16FrequencySeries(
                            "htilde: FD waveform", 
                            (LIGOTimeGPS){ 0, 0 }, 
                            0.0, 
                            deltaF, 
                            &lalStrainUnit, 
                            waveform_axes.strain.num_samples
                        );
                 
                    *hctilde = 
                        XLALCreateCOMPLEX16FrequencySeries(
                            "ctilde: FD waveform", 
                            (LIGOTimeGPS){ 0, 0 }, 
                            0.0, 
                            deltaF, 
                            &lalStrainUnit, 
                            waveform_axes.strain.num_samples
                        );
                }
                else // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
                {
                    *hptilde = 
                        XLALCreateCOMPLEX16FrequencySeries(
                            "htilde: FD waveform", 
                            (LIGOTimeGPS){ 0, 0 },
                            f_min, 
                            deltaF, 
                            &lalStrainUnit, 
                            waveform_axes.strain.num_samples
                        );
                    *hctilde = 
                        XLALCreateCOMPLEX16FrequencySeries(
                            "ctilde: FD waveform", 
                            (LIGOTimeGPS){ 0, 0 }, 
                            f_min, 
                            deltaF, 
                            &lalStrainUnit, 
                            waveform_axes.strain.num_samples
                        );
                }
             
                complex_strain_element_c *strain_fd = NULL;
                cudaToHost(
                    waveform_axes.strain.values, 
                    sizeof(complex_strain_element_c),
                    waveform_axes.strain.num_samples,
                    &strain_fd
                );


            for(j = 0; j < (*hptilde)->data->length; j++) 
            {    
                 (*hctilde)->data->data[j] = (complex double)strain_fd[j].cross;
                 (*hptilde)->data->data[j] = (complex double)strain_fd[j].plus;
            }
            break;
             
        case IMRPhenomXPHM:
             if(f_ref==0.0)
             {
                     /* Default reference frequency is minimum frequency */
                     f_ref = f_min;
             }

             /* Call the main waveform driver. Note that we pass the full spin vectors
                     with XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame being
                     effectively called in the initialization of the pPrec struct
            */

             ret = XLALSimIMRPhenomXPHM(
                     hptilde, hctilde,
                     m1, m2,
                     S1x, S1y, S1z,
                     S2x, S2y, S2z,
                     distance, inclination,
                     phiRef, f_min, f_max, deltaF, f_ref, LALparams
             );
                
             if (ret == XLAL_FAILURE)
             {
                     //XLAL_ERROR(XLAL_EFUNC);
             }

             break;

         default:
             XLALPrintError("FD version of approximant not implemented in lalsimulation\n");
             //XLAL_ERROR(XLAL_EINVAL);
         break;
     }
  
     double polariz=longAscNodes;
     if (polariz) {
       COMPLEX16 tmpP,tmpC;
       for (UINT4 idx=0;idx<(*hptilde)->data->length;idx++) 
       {
         tmpP=(*hptilde)->data->data[idx];
         tmpC=(*hctilde)->data->data[idx];
         (*hptilde)->data->data[idx] =cos(2.*polariz)*tmpP+sin(2.*polariz)*tmpC;
         (*hctilde)->data->data[idx]=cos(2.*polariz)*tmpC-sin(2.*polariz)*tmpP;
       }
     }
 
     //if (ret == XLAL_FAILURE) //XLAL_ERROR(XLAL_EFUNC);
     if (XLALSimInspiralWaveformParamsLookupEnableLIV(LALparams))
       ret = XLALSimLorentzInvarianceViolationTerm(hptilde, hctilde, m1/LAL_MSUN_SI, m2/LAL_MSUN_SI, distance, LALparams);
     //if (ret == XLAL_FAILURE) //XLAL_ERROR(XLAL_EFUNC);
  
     return ret;
}

int cuInspiralFD(
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
    companion_s companion_1,
    companion_s companion_2,
    double distance,                         /**< distance of source (m) */
    double inclination,                      /**< inclination of source (rad) */
    double phiRef,                           /**< reference orbital phase (rad) */
    double longAscNodes,                     /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    double eccentricity,                     /**< eccentricity at reference epoch */
    double meanPerAno,                       /**< mean anomaly of periastron */
    double deltaF,                           /**< sampling interval (Hz) */
    double f_min,                            /**< starting GW frequency (Hz) */
    double f_max,                            /**< ending GW frequency (Hz) */
    double f_ref,                            /**< Reference frequency (Hz) */
    LALDict *LALparams,                     /**< LAL dictionary containing accessory parameters */
    Approximant approximant                 /**< post-Newtonian approximant to use for waveform production */
    ) {  
    
    double m1  = companion_1.mass.kilograms; 
    double m2  = companion_2.mass.kilograms; 
    double S1x = companion_1.spin.x; 
    double S1y = companion_1.spin.y; 
    double S1z = companion_1.spin.z; 
    double S2x = companion_2.spin.x; 
    double S2y = companion_2.spin.y;
    double S2z = companion_2.spin.z;
    
    double chirplen, deltaT, f_nyquist;
    int chirplen_exp;
    int retval;
     size_t n;

    /* Apply condition that f_max rounds to the next power-of-two multiple
    * of deltaF.
    * Round f_max / deltaF to next power of two.
    * Set f_max to the new Nyquist frequency.
    * The length of the chirp signal is then 2 * f_nyquist / deltaF.
    * The time spacing is 1 / (2 * f_nyquist) */
    f_nyquist = f_max;
     if (deltaF != 0) {
         n = round(f_max / deltaF);
         if ((n & (n - 1))) { /* not a power of 2 */
                     frexp(n, &chirplen_exp);
                     f_nyquist = ldexp(1.0, chirplen_exp) * deltaF;
             XLAL_PRINT_WARNING("f_max/deltaF = %g/%g = %g is not a power of two: changing f_max to %g", f_max, deltaF, f_max/deltaF, f_nyquist);
     }
     }
    deltaT = 0.5 / f_nyquist;

    /* generate a FD waveform and condition it by applying tapers at
    * frequencies between a frequency below the requested f_min and
    * f_min; also wind the waveform in phase in case it would wrap-
    * around at the merger time */

    double tchirp, tmerge, textra, tshift;
    double fstart, fisco;
    double s;
    size_t k, k0, k1;

    /* if the requested low frequency is below the lowest Kerr ISCO
    * frequency then change it to that frequency */
    fisco = 1.0 / (pow(9.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    if (f_min > fisco)
     f_min = fisco;

    /* upper bound on the chirp time starting at f_min */
    tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, S1z, S2z);

    /* IMR model: estimate plunge and merger time */
    /* sometimes these waveforms have phases that
    * cause them to wrap-around an amount equal to
    * the merger-ringodwn time, so we will undo
    * that here */
    s = XLALSimInspiralFinalBlackHoleSpinBound(S1z, S2z);
    tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);

    /* new lower frequency to start the waveform: add some extra early
    * part over which tapers may be applied, the extra amount being
    * a fixed fraction of the chirp time; add some additional padding
    * equal to a few extra cycles at the low frequency as well for
    * safety and for other routines to use */
    textra = EXTRA_CYCLES / f_min;
    fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + EXTRA_TIME_FRACTION) * tchirp, m1, m2);

    /* revise (over-)estimate of chirp from new start frequency */
    tchirp = XLALSimInspiralChirpTimeBound(fstart, m1, m2, S1z, S2z);

    /* need a long enough segment to hold a whole chirp with some padding */
    /* length of the chirp in samples */
    chirplen = round((tchirp + tmerge + 2.0 * textra) / deltaT);
    /* make chirplen next power of two */
    frexp(chirplen, &chirplen_exp);
    chirplen = ldexp(1.0, chirplen_exp);
    /* frequency resolution */
    if (deltaF == 0.0)
     deltaF = 1.0 / (chirplen * deltaT);
    else if (deltaF > 1.0 / (chirplen * deltaT))
     XLAL_PRINT_WARNING("Specified frequency interval of %g Hz is too large for a chirp of duration %g s", deltaF, chirplen * deltaT);
    
    /* generate the waveform in the frequency domain starting at fstart */
    retval = cuInspiralChooseFDWaveform(hptilde, hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaF, fstart, f_max, f_ref, LALparams, approximant);
    
    //if (retval < 0)
     //XLAL_ERROR(XLAL_EFUNC);

    /* taper frequencies between fstart and f_min */
    k0 = round(fstart / (*hptilde)->deltaF);
    k1 = round(f_min / (*hptilde)->deltaF);
    /* make sure it is zero below fstart */
    for (k = 0; k < k0; ++k) {
     (*hptilde)->data->data[k] = 0.0;
     (*hctilde)->data->data[k] = 0.0;
    }
    /* taper between fstart and f_min */
    for ( ; k < k1; ++k) {
     double w = 0.5 - 0.5 * cos(M_PI * (k - k0) / (double)(k1 - k0));
     (*hptilde)->data->data[k] *= w;
     (*hctilde)->data->data[k] *= w;
    }
    /* make sure Nyquist frequency is zero */
    (*hptilde)->data->data[(*hptilde)->data->length - 1] = 0.0;
    (*hctilde)->data->data[(*hctilde)->data->length - 1] = 0.0;

    /* we want to make sure that this waveform will give something
    * sensible if it is later transformed into the time domain:
    * to avoid the end of the waveform wrapping around to the beginning,
    * we shift waveform backwards in time and compensate for this
    * shift by adjusting the epoch */
    tshift = round(tmerge / deltaT) * deltaT; /* integer number of time samples */
    for (k = 0; k < (*hptilde)->data->length; ++k) {
     double complex phasefac = cexp(2.0 * M_PI * I * k * deltaF * tshift);
     (*hptilde)->data->data[k] *= phasefac;
     (*hctilde)->data->data[k] *= phasefac;
    }
    XLALGPSAdd(&(*hptilde)->epoch, tshift);
    XLALGPSAdd(&(*hctilde)->epoch, tshift);

    return 0;
}

int cuInspiralTDFromFD(
     REAL8TimeSeries       **hplus,                    /**< +-polarization waveform */
     REAL8TimeSeries       **hcross,                   /**< x-polarization waveform */
     system_properties_s     system_properties,
     temporal_properties_s   temporal_properties,
     LALDict *LALparams,                         /**< LAL dictionary containing accessory parameters */
     Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
) {
    
    COMPLEX16FrequencySeries *hptilde = NULL;
    COMPLEX16FrequencySeries *hctilde = NULL;
    
    size_t end, k;
    
    // Unpack companion structs for readability:
    const companion_s companion_1 = system_properties.companion[0];
    const companion_s companion_2 = system_properties.companion[1];
    
    // Generate the conditioned waveform in the frequency domain note: redshift 
    // factor has already been applied above set deltaF = 0 to get a small
    // enough resolution:
    
   // printf("Cuda 3 | fstart: %f, f_ref %f, deltaT %f | \n", 
   //  temporal_properties.starting_frequency.hertz, 
   //  temporal_properties.reference_frequency.hertz, 
   //  temporal_properties.sampling_interval.seconds);
    int32_t return_value = 
        cuInspiralFD(
            &hptilde, 
            &hctilde, 
            system_properties.companion[0],
            system_properties.companion[1],
            system_properties.distance.meters,
            system_properties.inclination.radians, 
            system_properties.reference_orbital_phase.radians, 
            system_properties.ascending_node_longitude,
            system_properties.eccentricity, 
            system_properties.mean_periastron_anomaly,
            0.0, 
            temporal_properties.starting_frequency.hertz, 
            temporal_properties.ending_frequency.hertz, 
            temporal_properties.reference_frequency.hertz,
            LALparams,
            approximant
        );
    
    // we want to make sure that this waveform will give something
    // sensible if it is later transformed into the time domain:
    // to avoid the end of the waveform wrapping around to the beginning,
    // we shift waveform backwards in time and compensate for this
    // shift by adjusting the epoch -- note that XLALSimInspiralFD
    // guarantees that there is extra padding to do this 
     double tshift = round(temporal_properties.extra_time.seconds / temporal_properties.sampling_interval.seconds) * temporal_properties.sampling_interval.seconds; // integer number of samples 
     for (k = 0; k < hptilde->data->length; ++k) {
         double complex phasefac = cexp(2.0 * M_PI * I * k * hptilde->deltaF * tshift);
         hptilde->data->data[k] *= phasefac;
         hctilde->data->data[k] *= phasefac;
     }

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
                (1.0 + EXTRA_TIME_FRACTION)
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
        round(
            (new_inspiral_time_upper_bound.seconds 
            + temporal_properties.chirp_time_upper_bound.seconds)
            / temporal_properties.sampling_interval.seconds
        );
    
    tshift = 
          round(temporal_properties.extra_time.seconds / temporal_properties.sampling_interval.seconds) 
        * temporal_properties.sampling_interval.seconds; // Integer number of samples
    
    /* amount to snip off at the end is tshift */
    end = (*hplus)->data->length - round(tshift / temporal_properties.sampling_interval.seconds);

    /* snip off extra time at beginning and at the end */
    XLALResizeREAL8TimeSeries(*hplus, end - chirplen, chirplen);
    XLALResizeREAL8TimeSeries(*hcross, end - chirplen, chirplen);

    /* clean up */
    XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    XLALDestroyCOMPLEX16FrequencySeries(hctilde);
     
    return return_value;
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
    const frequencyUnit_t original_starting_frequency = starting_frequency; 
    
    // Init property structures:
    system_properties_s system_properties =
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
    
    temporal_properties_s temporal_properties =
        initTemporalProperties(
            sampling_interval,   // <-- Sampling interval (timeUnit_t).
            starting_frequency,  // <-- Starting GW frequency (frequencyUnit_t).
            reference_frequency, // <-- Reference GW frequency (frequencyUnit_t).
            system_properties,
            approximant
        );
    
    /* SEOBNR flag for spin aligned model version. 1 for SEOBNRv1, 2 for SEOBNRv2 */
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
    
    //printArrayE("First 10 2 Cuda", (*hplus)->data->data, 100);

    
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
                (*hplus)->data->length
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
        (*hplus)->data->length
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
          float64_2_t     **ret_strain
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
    double          ascending_node_longitude = 0.0;
    double          eccentricity             = 0.0;
    double          mean_periastron_anomaly  = 0.0;
    timeUnit_t      sampling_interval        = 
        initTimeSeconds(1.0/sample_rate.hertz);
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
    
    size_t new_array_size = (size_t)num_samples * sizeof(float64_2_t);

    float64_2_t *strain = (float64_2_t*)malloc(new_array_size);
    int32_t new_waveform_index = 0;
        for (int32_t index = 0; index < num_samples; index++) 
    {    
        new_waveform_index = waveform_num_samples - num_samples - 1 + index;
        strain[index].x = (float64_t)hplus->data->data[new_waveform_index];
        strain[index].y = (float64_t)hcross->data->data[new_waveform_index];
        }
    
    free(hcross->data->data); free(hplus->data->data);

    *ret_strain = strain;
}

#endif