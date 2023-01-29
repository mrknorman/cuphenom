#ifndef ZOMBIE_PHENOM_H
#define ZOMBIE_PHENOM_H

#include <omp.h>

#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#define LAL_PI 3.141592653589793238462643383279502884
#define LAL_MTSUN_SI 4.925491025543575903411922162094833998e-6
#define LAL_MSUN_SI 1.988409902147041637325262574352366540e30
#define LAL_MRSUN_SI 1.476625061404649406193430731479084713e3
#define LAL_PI_4 0.785398163397448309615660845819875721

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

void Nudge(
    double *x, 
    const double X, 
    const double epsilon
    ) {
    if (X != 0.0)
    {
        if (!gsl_fcmp(*x, X, epsilon))
        {
            XLAL_PRINT_INFO("Nudging value %.15g to %.15g\n", *x, X);
            *x = X;
        }
    }
    else
    {
        if (fabs(*x - X) < epsilon)
         *x = X;
    }
}

int ZombiePhenomDGenerateFD(
     COMPLEX16FrequencySeries **htilde,   // [out] FD waveform.
     const REAL8Sequence *freqs_in,       // Frequency points at which to evaluate the waveform (Hz)
     double deltaF,                       // If deltaF > 0, the frequency points given in freqs are uniformly spaced with
                                          // spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
                                          // Then we will use deltaF = 0 to create the frequency series we return.
     const double phi0,                   // phase at fRef
     const double fRef,                   // reference frequency [Hz]
     const double m1_in,                  // mass of companion 1 [solar masses]
     const double m2_in,                  // mass of companion 2 [solar masses]
     const double chi1_in,                // aligned-spin of companion 1
     const double chi2_in,                // aligned-spin of companion 2
     const double distance,               // distance to source (m)
     LALDict *extraParams,                // linked list containing the extra testing GR parameters
     NRTidal_version_type NRTidal_version // NRTidal version; either NRTidal_V or NRTidalv2_V or NoNRT_V in case of BBH baseline
) {
    timer_s timer;
    start_timer("Outside", &timer);
    
    timeUnit_t gps_time = initTimeSeconds(0.0f);
       
   // Make a pointer to LALDict to circumvent a memory leak
   // At the end we will check if we created a LALDict in extraParams
   // and destroy it if we did:
   
   LALDict *extraParams_in = extraParams;
   REAL8Sequence *amp_tidal = NULL; /* Tidal amplitude series; required only for IMRPhenomD_NRTidalv2 */
   double dquadmon1_in = 0., dquadmon2_in = 0., lambda1_in = 0, lambda2_in = 0.;
   if (NRTidal_version == NRTidalv2_V) 
   {
     dquadmon1_in = XLALSimInspiralWaveformParamsLookupdQuadMon1(extraParams);
     dquadmon2_in = XLALSimInspiralWaveformParamsLookupdQuadMon2(extraParams);
     lambda1_in = XLALSimInspiralWaveformParamsLookupTidalLambda1(extraParams);
     lambda2_in = XLALSimInspiralWaveformParamsLookupTidalLambda2(extraParams);
   }
  
   double chi1, chi2, m1, m2, dquadmon1, dquadmon2, lambda1, lambda2;
   if (m1_in>=m2_in) 
   {
        chi1      = chi1_in;
        chi2      = chi2_in;
        m1        = m1_in;
        m2        = m2_in;
        dquadmon1 = dquadmon1_in;
        dquadmon2 = dquadmon2_in;
        lambda1   = lambda1_in;
        lambda2   = lambda2_in;
    } 
    else  // swap spins and masses
    {
        chi1      = chi2_in;
        chi2      = chi1_in;
        m1        = m2_in;
        m2        = m1_in;
        dquadmon1 = dquadmon2_in;
        dquadmon2 = dquadmon1_in;
        lambda1   = lambda2_in;
        lambda2   = lambda1_in;
        if (NRTidal_version == NRTidalv2_V) 
        {
            XLALSimInspiralWaveformParamsInsertdQuadMon1(extraParams, dquadmon1);
            XLALSimInspiralWaveformParamsInsertdQuadMon2(extraParams, dquadmon2);
        }
    }
     
    print_timer("1", &timer);

    int status = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initiate useful powers of pi.");
  
   /* Find frequency bounds */
   if (!freqs_in || !freqs_in->data) XLAL_ERROR(XLAL_EFAULT);
   double f_min = freqs_in->data[0];
   double f_max = freqs_in->data[freqs_in->length - 1];
   XLAL_CHECK(f_min > 0, XLAL_EDOM, "Minimum frequency must be positive.\n");
   XLAL_CHECK(f_max >= 0, XLAL_EDOM, "Maximum frequency must be non-negative.\n");
  
   const double M = m1 + m2;
   double eta = m1 * m2 / (M * M);
  
   if (eta > 0.25)
       Nudge(&eta, 0.25, 1e-6);
   if (eta > 0.25 || eta < 0.0)
       XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");
  
   const double M_sec = M * LAL_MTSUN_SI;
  
   /* Compute the amplitude pre-factor */
   const double amp0 = 2. * sqrt(5. / (64.*LAL_PI)) * M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;
    
     
   print_timer("2", &timer);

   size_t num_frequency_samples = 0;
   UINT4 offset = 0; // Index shift between freqs and the frequency series
   REAL8Sequence *freqs = NULL;
     
    // Freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0:
    if (deltaF > 0)  
    { 
        // Set up output array with size closest power of 2:
        num_frequency_samples = calcNextPow2(f_max / deltaF) + 1;
       
        // Coalesce at gps_time = 0:
        // Shift by overall length in time:  
        gps_time = addTimes(2, gps_time, initTimeSeconds(-1. / deltaF));
        
         *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", (LIGOTimeGPS){ 0, 0 }, 0.0, deltaF, &lalStrainUnit, num_frequency_samples);

        // Recreate freqs using only the lower and upper bounds
        size_t iStart = (size_t) (f_min / deltaF);
        size_t iStop = (size_t) (f_max / deltaF);
       
     XLAL_CHECK ( (iStop<=num_frequency_samples) && (iStart<=iStop), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, num_frequency_samples);
     freqs = XLALCreateREAL8Sequence(iStop - iStart);
     if (!freqs)
       XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
     for (size_t i = iStart; i < iStop; i++)
       freqs->data[i-iStart] = i*deltaF;
     offset = iStart;
   } else { // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
     num_frequency_samples = freqs_in->length;
     *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", (LIGOTimeGPS){ 0, 0 }, f_min, deltaF, &lalStrainUnit, num_frequency_samples);
     XLAL_CHECK ( *htilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu from sequence.", num_frequency_samples);
     offset = 0;
     freqs = XLALCreateREAL8Sequence(freqs_in->length);
     if (!freqs)
       XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
     for (size_t i=0; i<freqs_in->length; i++)
       freqs->data[i] = freqs_in->data[i];
   }
     
    print_timer("3", &timer);
  
   memset((*htilde)->data->data, 0, num_frequency_samples * sizeof(COMPLEX16));
    
   XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);
  
   // Calculate phenomenological parameters
   const double finspin = ZombFinalSpin0815(eta, chi1, chi2); //FinalSpin0815 - 0815 is like a version number
  
   if (finspin < MIN_FINAL_SPIN)
           XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system_properties are small, \
                          the model might misbehave here.", finspin);
  
   IMRPhenomDAmplitudeCoefficients *pAmp;
   pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
   ZombComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1, chi2, finspin);
   if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
   if (extraParams==NULL)
     extraParams=XLALCreateDict();
   XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
   IMRPhenomDPhaseCoefficients *pPhi;
   pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
   ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1, chi2, finspin, extraParams);
    
   if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
   PNPhasingSeries *pn = NULL;
   XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1, chi2, extraParams);
   if (!pn) XLAL_ERROR(XLAL_EFUNC);
  
   // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
   // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
   double testGRcor=1.0;
   testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);
  
   // was not available when PhenomD was tuned.
   pn->v[6] -= (Subtract3PNSS(m1, m2, M, eta, chi1, chi2) * pn->v[0]) * testGRcor;
  
   PhiInsPrefactors phi_prefactors;
   status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
   XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");
     
    print_timer("4", &timer);
  
   // Compute coefficients to make phase C^1 continuous (phase and first derivative)
   ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, 1.0, 1.0);
  
   //time shift so that peak amplitude is approximately at t=0
   //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
   const double t0 = DPhiMRD(pAmp->fmaxCalc, pPhi, 1.0, 1.0);
  
   AmpInsPrefactors amp_prefactors;
   status = init_amp_ins_prefactors(&amp_prefactors, pAmp);
   XLAL_CHECK(XLAL_SUCCESS == status, status, "init_amp_ins_prefactors failed");
  
   // incorporating fRef
   const double MfRef = M_sec * fRef;
   UsefulPowers powers_of_fRef;
   status = init_useful_powers(&powers_of_fRef, MfRef);
   XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers failed for MfRef");
    
   const double phifRef = IMRPhenDPhase(MfRef, pPhi, pn, &powers_of_fRef, &phi_prefactors, 1.0, 1.0);
  
   // factor of 2 b/c phi0 is orbital phase
   const double phi_precalc = 2.*phi0 + phifRef;
  
   int status_in_for = XLAL_SUCCESS;
   int ret = XLAL_SUCCESS;
   /* Now generate the waveform */
     
    print_timer("5", &timer);

    
   if (NRTidal_version == NRTidalv2_V) {
     /* Generate the tidal amplitude (Eq. 24 of arxiv: 1905.06011) to add to BBH baseline; only for IMRPhenomD_NRTidalv2 */
     amp_tidal = XLALCreateREAL8Sequence(freqs->length);
     ret = XLALSimNRTunedTidesFDTidalAmplitudeFrequencySeries(amp_tidal, freqs, m1, m2, lambda1, lambda2);
     XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to generate tidal amplitude series to construct IMRPhenomD_NRTidalv2 waveform.");
     /* Generated tidal amplitude corrections */
    #pragma omp parallel for
     for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
       double Mf = M_sec * freqs->data[i];
       double ampT = amp_tidal->data[i];
       int j = i + offset; // shift index for frequency series if needed
  
       UsefulPowers powers_of_f;
       status_in_for = init_useful_powers(&powers_of_f, Mf);
       if (XLAL_SUCCESS != status_in_for)
       {
         XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
         status = status_in_for;
       }
       else {
         double amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
        double phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors, 1.0, 1.0);
         phi -= t0*(Mf-MfRef) + phi_precalc;
         ((*htilde)->data->data)[j] = amp0 * (amp+2*sqrt(LAL_PI/5.)*ampT) * cexp(-I * phi);
       }
     }
   } 
    
    else {
      #pragma omp parallel for
        
       for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
       
       double Mf = M_sec * freqs->data[i];
       int j = i + offset; // shift index for frequency series if needed
  
       UsefulPowers powers_of_f;
       status_in_for = init_useful_powers(&powers_of_f, Mf);
       if (XLAL_SUCCESS != status_in_for)
       {
         XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
         status = status_in_for;
       }
       else {
         double amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
         double phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors, 1.0, 1.0);
         phi -= t0*(Mf-MfRef) + phi_precalc;
         ((*htilde)->data->data)[j] = amp0 * amp * cexp(-I * phi); 
       }
     }
   }
     
   print_timer("6", &timer);

  
   LALFree(pAmp);
   LALFree(pPhi);
   LALFree(pn);
   XLALDestroyREAL8Sequence(freqs);
   XLALDestroyREAL8Sequence(amp_tidal);
  
   /* If extraParams was allocated in this function and not passed in
   * we need to free it to prevent a leak */
   if (extraParams && !extraParams_in) {
     XLALDestroyDict(extraParams);
   } else {
     XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_ALL);
   }
  
   return status;
}

int cuPhenomDGenerateFD(
     COMPLEX16FrequencySeries **htilde,   // [out] FD waveform.
     const REAL8Sequence *freqs_in,       // Frequency points at which to evaluate the waveform (Hz)
     double deltaF,                       // If deltaF > 0, the frequency points given in freqs are uniformly spaced with
                                          // spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
                                          // Then we will use deltaF = 0 to create the frequency series we return.
     const double phi0,                   // phase at fRef
     const double fRef,                   // reference frequency [Hz]
     const double m1_in,                  // mass of companion 1 [solar masses]
     const double m2_in,                  // mass of companion 2 [solar masses]
     const double chi1_in,                // aligned-spin of companion 1
     const double chi2_in,                // aligned-spin of companion 2
     const double distance,               // distance to source (m)
     LALDict *extraParams,                // linked list containing the extra testing GR parameters
     NRTidal_version_type NRTidal_version // NRTidal version; either NRTidal_V or NRTidalv2_V or NoNRT_V in case of BBH baseline
) {
    timer_s timer;
    start_timer("Outside", &timer);

    timeUnit_t gps_time = initTimeSeconds(0.0f);

    // Make a pointer to LALDict to circumvent a memory leak
    // At the end we will check if we created a LALDict in extraParams
    // and destroy it if we did:

    LALDict *extraParams_in = extraParams;
    REAL8Sequence *amp_tidal = NULL; /* Tidal amplitude series; required only for IMRPhenomD_NRTidalv2 */
    double dquadmon1_in = 0., dquadmon2_in = 0., lambda1_in = 0, lambda2_in = 0.;
    if (NRTidal_version == NRTidalv2_V) 
    {
        dquadmon1_in = XLALSimInspiralWaveformParamsLookupdQuadMon1(extraParams);
        dquadmon2_in = XLALSimInspiralWaveformParamsLookupdQuadMon2(extraParams);
        lambda1_in = XLALSimInspiralWaveformParamsLookupTidalLambda1(extraParams);
        lambda2_in = XLALSimInspiralWaveformParamsLookupTidalLambda2(extraParams);
    }

    double chi1, chi2, m1, m2, dquadmon1, dquadmon2, lambda1, lambda2;
    if (m1_in>=m2_in) 
    {
        chi1      = chi1_in;
        chi2      = chi2_in;
        m1        = m1_in;
        m2        = m2_in;
        dquadmon1 = dquadmon1_in;
        dquadmon2 = dquadmon2_in;
        lambda1   = lambda1_in;
        lambda2   = lambda2_in;
    } 
    else  // swap spins and masses
    {
        chi1      = chi2_in;
        chi2      = chi1_in;
        m1        = m2_in;
        m2        = m1_in;
        dquadmon1 = dquadmon2_in;
        dquadmon2 = dquadmon1_in;
        lambda1   = lambda2_in;
        lambda2   = lambda1_in;
        if (NRTidal_version == NRTidalv2_V) 
        {
            XLALSimInspiralWaveformParamsInsertdQuadMon1(extraParams, dquadmon1);
            XLALSimInspiralWaveformParamsInsertdQuadMon2(extraParams, dquadmon2);
        }
    }

    print_timer("1", &timer);

    int status = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initiate useful powers of pi.");

    /* Find frequency bounds */
    if (!freqs_in || !freqs_in->data) XLAL_ERROR(XLAL_EFAULT);
    
    double f_min = freqs_in->data[0];
    double f_max = freqs_in->data[freqs_in->length - 1];
    XLAL_CHECK(f_min > 0, XLAL_EDOM, "Minimum frequency must be positive.\n");
    XLAL_CHECK(f_max >= 0, XLAL_EDOM, "Maximum frequency must be non-negative.\n");

    const double M = m1 + m2;
    double eta = m1 * m2 / (M * M);

    if (eta > 0.25)
        Nudge(&eta, 0.25, 1e-6);
    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    const double M_sec = M * LAL_MTSUN_SI;

    /* Compute the amplitude pre-factor */
    const double amp0 = 2. * sqrt(5. / (64.*LAL_PI)) * M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;


    print_timer("2", &timer);

    size_t num_frequency_samples = 0;
    int32_t offset = 0; // Index shift between freqs and the frequency series
    REAL8Sequence *freqs = NULL;
    
    
    
    
    
    
    if (deltaF > 0) // Freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0:
    { 
        // Set up output array with size closest power of 2:
        num_frequency_samples = calcNextPow2(f_max / deltaF) + 1;
        *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", (LIGOTimeGPS){ 0, 0 }, 0.0, deltaF, &lalStrainUnit, num_frequency_samples);
        
        // Coalesce at gps_time = 0:
        // Shift by overall length in time:  
        gps_time = addTimes(2, gps_time, initTimeSeconds(-1. / deltaF));
    }
    else // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
    {
        num_frequency_samples = freqs_in->length;
        *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", (LIGOTimeGPS){ 0, 0 }, f_min, deltaF, &lalStrainUnit, num_frequency_samples);
        freqs->data = freqs_in->data;
        offset = 0;
    }
    
    complex float *strain_fd      = NULL;
    complex float *frequency_axis = NULL;
    
    if (deltaF > 0)  
    { 




        // Recreate freqs using only the lower and upper bounds
        size_t iStart = (size_t) (f_min / deltaF);
        size_t iStop = (size_t) (f_max / deltaF);

        XLAL_CHECK ( (iStop<=num_frequency_samples) && (iStart<=iStop), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, num_frequency_samples);
        freqs = XLALCreateREAL8Sequence(iStop - iStart);
        for (size_t i = iStart; i < iStop; i++)
            freqs->data[i-iStart] = i*deltaF;
        offset = iStart;
    } 

    print_timer("3", &timer);

    memset((*htilde)->data->data, 0, num_frequency_samples * sizeof(COMPLEX16));

    XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);

    // Calculate phenomenological parameters
    const double finspin = ZombFinalSpin0815(eta, chi1, chi2); //FinalSpin0815 - 0815 is like a version number

    if (finspin < MIN_FINAL_SPIN)
    XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system_properties are small, \
          the model might misbehave here.", finspin);

    IMRPhenomDAmplitudeCoefficients *pAmp;
    pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
    ZombComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1, chi2, finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    if (extraParams==NULL)
    extraParams=XLALCreateDict();
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1, chi2, finspin, extraParams);

    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1, chi2, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    double testGRcor=1.0;
    testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);

    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSS(m1, m2, M, eta, chi1, chi2) * pn->v[0]) * testGRcor;

    PhiInsPrefactors phi_prefactors;
    status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

    print_timer("4", &timer);

    // Compute coefficients to make phase C^1 continuous (phase and first derivative)
    ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, 1.0, 1.0);

    //time shift so that peak amplitude is approximately at t=0
    //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
    const double t0 = DPhiMRD(pAmp->fmaxCalc, pPhi, 1.0, 1.0);

    AmpInsPrefactors amp_prefactors;
    status = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_amp_ins_prefactors failed");

    // incorporating fRef
    const double MfRef = M_sec * fRef;
    UsefulPowers powers_of_fRef;
    status = init_useful_powers(&powers_of_fRef, MfRef);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers failed for MfRef");

    const double phifRef = IMRPhenDPhase(MfRef, pPhi, pn, &powers_of_fRef, &phi_prefactors, 1.0, 1.0);

    // factor of 2 b/c phi0 is orbital phase
    const double phi_precalc = 2.*phi0 + phifRef;

    int status_in_for = XLAL_SUCCESS;
    int ret = XLAL_SUCCESS;
    /* Now generate the waveform */

    print_timer("5", &timer);

    if (NRTidal_version == NRTidalv2_V) 
    {
        /* Generate the tidal amplitude (Eq. 24 of arxiv: 1905.06011) to add to BBH baseline; only for IMRPhenomD_NRTidalv2 */
        amp_tidal = XLALCreateREAL8Sequence(freqs->length);
        ret = XLALSimNRTunedTidesFDTidalAmplitudeFrequencySeries(amp_tidal, freqs, m1, m2, lambda1, lambda2);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to generate tidal amplitude series to construct IMRPhenomD_NRTidalv2 waveform.");
        /* Generated tidal amplitude corrections */
        #pragma omp parallel for
        for (UINT4 i=0; i<freqs->length; i++) 
        { // loop over frequency points in sequence
            double Mf = M_sec * freqs->data[i];
            double ampT = amp_tidal->data[i];
            int j = i + offset; // shift index for frequency series if needed

            UsefulPowers powers_of_f;
            status_in_for = init_useful_powers(&powers_of_f, Mf);
            if (XLAL_SUCCESS != status_in_for)
            {
                XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
                status = status_in_for;
            }
            else 
            {
                double amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
                double phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors, 1.0, 1.0);
                phi -= t0*(Mf-MfRef) + phi_precalc;
                ((*htilde)->data->data)[j] = amp0 * (amp+2*sqrt(LAL_PI/5.)*ampT) * cexp(-I * phi);
            }
        }
    } 
    else {
        #pragma omp parallel for
        for (UINT4 i=0; i<freqs->length; i++) 
        { // loop over frequency points in sequence
            double Mf = M_sec * freqs->data[i];
            int j = i + offset; // shift index for frequency series if needed

            UsefulPowers powers_of_f;
            status_in_for = init_useful_powers(&powers_of_f, Mf);
            if (XLAL_SUCCESS != status_in_for)
            {
                XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
                status = status_in_for;
            }
            else 
            {
                double amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
                double phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors, 1.0, 1.0);
                phi -= t0*(Mf-MfRef) + phi_precalc;
                ((*htilde)->data->data)[j] = amp0 * amp * cexp(-I * phi); 
            }
        }
    }

    print_timer("6", &timer);


    LALFree(pAmp);
    LALFree(pPhi);
    LALFree(pn);
    XLALDestroyREAL8Sequence(freqs);
    XLALDestroyREAL8Sequence(amp_tidal);
    
    return status;
}

int ZombieIMRPhenomDGenerateFD(
     COMPLEX16FrequencySeries **htilde, /**< [out] FD waveform */
     const double phi0,                  /**< Orbital phase at fRef (rad) */
     const double fRef_in,               /**< reference frequency (Hz) */
     const double deltaF,                /**< Sampling frequency (Hz) */
     const double m1_SI,                 /**< Mass of companion 1 (kg) */
     const double m2_SI,                 /**< Mass of companion 2 (kg) */
     const double chi1,                  /**< Aligned-spin parameter of companion 1 */
     const double chi2,                  /**< Aligned-spin parameter of companion 2 */
     const double f_min,                 /**< Starting GW frequency (Hz) */
     const double f_max,                 /**< End frequency; 0 defaults to Mf = \ref f_CUT */
     const double distance,               /**< Distance of source (m) */
     LALDict *extraParams, /**< linked list containing the extra testing GR parameters */
     NRTidal_version_type  NRTidal_version /**< Version of NRTides; can be one of NRTidal versions or NoNRT_V for the BBH baseline */
) {
    /* external: SI; internal: solar masses */
    const double m1 = m1_SI / LAL_MSUN_SI;
    const double m2 = m2_SI / LAL_MSUN_SI;

    /* check inputs for sanity */
    XLAL_CHECK(0 != htilde, XLAL_EFAULT, "htilde is null");
    if (*htilde) XLAL_ERROR(XLAL_EFAULT);
    if (fRef_in < 0) XLAL_ERROR(XLAL_EDOM, "fRef_in must be positive (or 0 for 'ignore')\n");
    if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
    if (m1 <= 0) XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
    if (m2 <= 0) XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");
    if (f_min <= 0) XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
    if (f_max < 0) XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");
    if (distance <= 0) XLAL_ERROR(XLAL_EDOM, "distance must be positive\n");

    const double q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

    if (q > 1000)
     XLAL_PRINT_WARNING("Warning: The model is not supported for high mass ratio, see MAX_ALLOWED_MASS_RATIO\n");

    if (chi1 > 1.0 || chi1 < -1.0 || chi2 > 1.0 || chi2 < -1.0)
     XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    // if no reference frequency given, set it to the starting GW frequency
    double fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    const double M_sec = (m1+m2) * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
    const double fCut = f_CUT/M_sec; // convert Mf -> Hz
    // Somewhat arbitrary end point for the waveform.
    // Chosen so that the end of the waveform is well after the ringdown.
    if (fCut <= f_min)
     XLAL_ERROR(XLAL_EDOM, "(fCut = %g Hz) <= f_min = %g\n", fCut, f_min);

     /* default f_max to Cut */
    double f_max_prime = f_max;
    f_max_prime = f_max ? f_max : fCut;
    f_max_prime = (f_max_prime > fCut) ? fCut : f_max_prime;
    if (f_max_prime <= f_min)
     XLAL_ERROR(XLAL_EDOM, "f_max <= f_min\n");

    // Use fLow, fHigh, deltaF to compute freqs sequence
    // Instead of building a full sequency we only transfer the boundaries and let
    // the internal core function do the rest (and properly take care of corner cases).
    REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
    freqs->data[0] = f_min;
    freqs->data[1] = f_max_prime;
    
    int status = 
        cuPhenomDGenerateFD(
            htilde, 
            freqs, 
            deltaF, 
            phi0, 
            fRef,
            m1, 
            m2, 
            chi1, 
            chi2,
            distance, 
            extraParams, 
            NRTidal_version
        );
    
    XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to generate IMRPhenomD waveform.");
    XLALDestroyREAL8Sequence(freqs);

    if (f_max_prime < f_max) {
     // The user has requested a higher f_max than Mf=fCut.
     // Resize the frequency series to fill with zeros beyond the cutoff frequency.
     size_t n = (*htilde)->data->length;
     size_t n_full = calcNextPow2(f_max / deltaF) + 1; // we actually want to have the length be a power of 2 + 1
     *htilde = XLALResizeCOMPLEX16FrequencySeries(*htilde, 0, n_full);
     XLAL_CHECK ( *htilde, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, fCut, n_full, f_max );
    }

    return XLAL_SUCCESS;
}

int ZombieSimInspiralChooseFDWaveform(
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
     )
{
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
  
     switch (approximant)
     {
         case IMRPhenomD:
             // Call the waveform driver routine:
             ret = ZombieIMRPhenomDGenerateFD(hptilde, phiRef, f_ref, deltaF, m1, m2,
                   S1z, S2z, f_min, f_max, distance, LALparams, 4);
             if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
             /* Produce both polarizations */
             *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                     &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                     &((*hptilde)->sampleUnits), (*hptilde)->data->length);
             for(j = 0; j < (*hptilde)->data->length; j++) {
                 (*hctilde)->data->data[j] = -I*cfac * (*hptilde)->data->data[j];
                 (*hptilde)->data->data[j] *= pfac;
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
                     XLAL_ERROR(XLAL_EFUNC);
             }

             break;

         default:
             XLALPrintError("FD version of approximant not implemented in lalsimulation\n");
             XLAL_ERROR(XLAL_EINVAL);
     }
  
     double polariz=longAscNodes;
     if (polariz) {
       COMPLEX16 tmpP,tmpC;
       for (UINT4 idx=0;idx<(*hptilde)->data->length;idx++) {
         tmpP=(*hptilde)->data->data[idx];
         tmpC=(*hctilde)->data->data[idx];
         (*hptilde)->data->data[idx] =cos(2.*polariz)*tmpP+sin(2.*polariz)*tmpC;
         (*hctilde)->data->data[idx]=cos(2.*polariz)*tmpC-sin(2.*polariz)*tmpP;
       }
     }
  
     if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
     if (XLALSimInspiralWaveformParamsLookupEnableLIV(LALparams))
       ret = XLALSimLorentzInvarianceViolationTerm(hptilde, hctilde, m1/LAL_MSUN_SI, m2/LAL_MSUN_SI, distance, LALparams);
     if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
  
     return ret;
}

int InspiralFD(
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
    retval = ZombieSimInspiralChooseFDWaveform(hptilde, hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaF, fstart, f_max, f_ref, LALparams, approximant);
    if (retval < 0)
     XLAL_ERROR(XLAL_EFUNC);

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

int InspiralTDFromFD(
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
    int32_t return_value = 
        InspiralFD(
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

    //XLALGPSAdd(&hptilde->epoch, tshift);
    //XLALGPSAdd(&hctilde->epoch, tshift);

    // Rransform the waveform into the time domain:
    size_t chirplen = 2 * (hptilde->data->length - 1);
    *hplus = XLALCreateREAL8TimeSeries("H_PLUS", &hptilde->epoch, 0.0, temporal_properties.sampling_interval.seconds, &lalStrainUnit, chirplen);
    *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &hctilde->epoch, 0.0, temporal_properties.sampling_interval.seconds, &lalStrainUnit, chirplen);
    
    size_t num_waveform_samples = 2 * (hptilde->data->length - 1);
        performIRFFT(
            hptilde->data->data,
            hctilde->data->data,
            (*hplus)->data->data,
            (*hcross)->data->data,
            temporal_properties,
            initFrequencyHertz(hptilde->deltaF),
            (int32_t)num_waveform_samples
        );

    /* compute how long a chirp we should have */
    /* revised estimate of chirp length from new start frequency */
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
    
    double tshift = 
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
    int32_t return_value = InspiralTDFromFD(
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
    XLALSimInspiralTDConditionStage1(
        *hplus, 
        *hcross, 
        EXTRA_TIME_FRACTION*temporal_properties.chirp_time_upper_bound.seconds
        + temporal_properties.extra_time.seconds, 
        original_starting_frequency.hertz
    );
    
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