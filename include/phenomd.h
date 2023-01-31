#include "phenomd_data.h"
#include "phenomd_structures.h"

static inline double Square(
    const double number
) {
    return number*number;
}

static inline double Cube(
    const double number)
{
    return number*number*number;
}

static inline double Quart(double number)
{
    double pow2 = Square(number);
    return pow2 * pow2;
}

inline size_t calcNextPow2(const double n)
{
   // Use pow here, not bit-wise shift, as the latter seems to run against an
   // upper cutoff long before SIZE_MAX, at least on some platforms:
   return (size_t) pow(2,ceil(log2(n)));
}

// See phenomd_data.h for the coefficient terms.
inline double calcCoefficient(
    const double *terms,
    const double  eta,
    const double  eta2, 
    const double  xi
    ) {
    
   return 
      terms[0] + terms[1]*eta
   + (terms[2] + terms[3]*eta + terms[ 4]*eta2
   + (terms[5] + terms[6]*eta + terms[ 7]*eta2)*xi
   + (terms[8] + terms[9]*eta + terms[10]*eta2)*xi*xi)*xi;
}

void *calcCoefficients(
    const double **terms,
    const size_t   num_terms,
    const double   eta,
    const double   eta2, 
    const double   xi,
          double  *coeffcients
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

/**
* Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
* (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
* was not available when PhenomD was tuned.
*/
double Subtract3PNSS(double m1, double m2, double M, double eta, double chi1, double chi2){
   long double m1M = (long double) (m1 / M);
   long double m2M = (long double) (m2 / M);
   double pn_ss3 = (double) ((326.75L/1.12L + 557.5L/1.8L*(long double)eta)*(long double)(eta*chi1*chi2));
   pn_ss3 += 
       (double)(((4703.5L/8.4L+2935.L/6.L*m1M-120.L*m1M*m1M) + (-4108.25L/6.72L-108.5L/1.2L*m1M+125.5L/3.6L*m1M*m1M))*m1M*m1M 
       *(long double)(chi1*chi1));
   pn_ss3 += 
       (double)(((4703.5L/8.4L+2935.L/6.L*m2M-120.L*m2M*m2M) + (-4108.25L/6.72L-108.5L/1.2L*m2M+125.5L/3.6L*m2M*m2M))*m2M*m2M 
       *(long double)(chi2*chi2));
   return pn_ss3;
}

double chiPN(double Seta, double eta, double chi1, double chi2)
{
   // Convention m1 >= m2 and chi1 is the spin on m1
   // The 0.5 factor missing in the definitions of chi_s and chi_a is
   // recovered in the return expresion
   double chi_s = (chi1 + chi2);
   double chi_a = (chi1 - chi2);
  
   return 0.5 * (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a);
}

static double ZombEradRational0815_s(double eta, double s)
{
    double eta2 = eta * eta;
    double eta3 = eta2 * eta;
  
     return (eta * (0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2 + 3.352411249771192 * eta3) *
             (1. + (-0.0030302335878845507 - 2.0066110851351073 * eta + 7.7050567802399215 * eta2) * s)) /
            (1. + (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * s);
}

double ZombPhenomInternal_EradRational0815(double eta, double chi1, double chi2)
{
     // Convention m1 >= m2
     double Seta = sqrt(1.0 - 4.0 * eta);
     double m1 = 0.5 * (1.0 + Seta);
     double m2 = 0.5 * (1.0 - Seta);
     double m1s = m1 * m1;
     double m2s = m2 * m2;
     // arXiv:1508.07250
     double s = (m1s * chi1 + m2s * chi2) / (m1s + m2s);
  
     return ZombEradRational0815_s(eta, s);
}

double Zombfring(double eta, double chi1, double chi2, double finspin) {
   double return_val;
  
   if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");
  
   gsl_interp_accel *acc = gsl_interp_accel_alloc();
   gsl_spline *iFring = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
   gsl_spline_init(iFring, QNMData_a, QNMData_fring, QNMData_length);
  
   return_val = gsl_spline_eval(iFring, finspin, acc) / (1.0 - ZombPhenomInternal_EradRational0815(eta, chi1, chi2));
  
   gsl_spline_free(iFring);
   gsl_interp_accel_free(acc);
    
   return return_val;
}

 double Zombfdamp(double eta, double chi1, double chi2, double finspin) {
   double return_val;
  
   if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fdamp function: final spin > 1.0 not supported\n");
  
   gsl_interp_accel *acc = gsl_interp_accel_alloc();
   gsl_spline *iFdamp = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
   gsl_spline_init(iFdamp, QNMData_a, QNMData_fdamp, QNMData_length);
  
   return_val = gsl_spline_eval(iFdamp, finspin, acc) / (1.0 - ZombPhenomInternal_EradRational0815(eta, chi1, chi2));
    
   gsl_spline_free(iFdamp);
   gsl_interp_accel_free(acc);
   return return_val;
}

typedef enum tagNRTidal_version_type {
  NRTidal_V, /**< version NRTidal: based on https://arxiv.org/pdf/1706.02969.pdf*/
  NRTidalv2_V, /**< version NRTidalv2: https://arxiv.org/abs/1905.06011 */
  NRTidalv2NoAmpCorr_V, /**< version NRTidalv2, without amplitude corrections */
  NRTidalv2NSBH_V, /**< version NRTidalv2: https://arxiv.org/abs/1905.06011 with amplitude corrections for NSBH (used for SEOBNRv4ROM_NRTidalv2_NSBH) */
  NoNRT_V /**< special case for PhenomPv2 BBH baseline */
} NRTidal_version_type;


useful_powers_s powers_of_pi;   

///////////////////////////// Amplitude: Intermediate functions ////////////////////////
 
// Phenom coefficients delta0, ..., delta4 determined from collocation method
// (constraining 3 values and 2 derivatives)
// AmpIntAnsatzFunc[]

// The Newtonian term in LAL is fine and we should use exactly the same (either hardcoded or call).
// We just use the Mathematica expression for convenience.
/**
 * Inspiral amplitude plus rho phenom coefficents. rho coefficients computed
 * in rho1_fun, rho2_fun, rho3_fun functions.
 * Amplitude is a re-expansion. See 1508.07253 and Equation 29, 30 and Appendix B arXiv:1508.07253 for details
 */
 double AmpInsAnsatz(double Mf, useful_powers_s * powers_of_Mf, AmpInsPrefactors * prefactors) {
  
   return 1 + powers_of_Mf->two_thirds * prefactors->two_thirds
       + powers_of_Mf->four_thirds * prefactors->four_thirds
                         + powers_of_Mf->five_thirds * prefactors->five_thirds
                         + powers_of_Mf->seven_thirds * prefactors->seven_thirds + powers_of_Mf->eight_thirds * prefactors->eight_thirds
                         + Mf * (prefactors->one + Mf * prefactors->two + powers_of_Mf->two * prefactors->three);
}

/**
 * Ansatz for the merger-ringdown amplitude. Equation 19 arXiv:1508.07253
 */
 double AmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
   double fRD = p->fRD;
   double fDM = p->fDM;
   double gamma1 = p->merger_ringdown[0];
   double gamma2 = p->merger_ringdown[1];
   double gamma3 = p->merger_ringdown[2];
   double fDMgamma3 = fDM*gamma3;
   double fminfRD = f - fRD;
   return exp( -(fminfRD)*gamma2 / (fDMgamma3) )
     * (fDMgamma3*gamma1) / (Square(fminfRD) + Square(fDMgamma3));
}
 
// Ansatz for the intermediate amplitude. Equation 21 arXiv:1508.07253:
double AmpIntAnsatz(
    const double                           Mf, 
    const IMRPhenomDAmplitudeCoefficients *amplitude_coefficients) {
    double Mf2 = Square(Mf);
    
    const double delta0 = amplitude_coefficients->intermediate[0];
    const double delta1 = amplitude_coefficients->intermediate[1];
    const double delta2 = amplitude_coefficients->intermediate[2];
    const double delta3 = amplitude_coefficients->intermediate[3];
    const double delta4 = amplitude_coefficients->intermediate[4];
       
    return delta0 + Mf*delta1 + Mf2*(delta2 + Mf*delta3 + Mf2*delta4);
}


// Call ComputeIMRPhenomDAmplitudeCoefficients() first!
/**
 * This function computes the IMR amplitude given phenom coefficients.
 * Defined in VIII. Full IMR Waveforms arXiv:1508.07253
 */
 double IMRPhenDAmplitude(double f, IMRPhenomDAmplitudeCoefficients *p, useful_powers_s *powers_of_f, AmpInsPrefactors * prefactors) {
   // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
   // The inspiral, intermediate and merger-ringdown amplitude parts
  
   // Transition frequencies
   p->fInsJoin = AMP_fJoin_INS;
   p->fMRDJoin = p->fmaxCalc;
  
   double AmpPreFac = prefactors->amp0 * powers_of_f->m_seven_sixths;
  
   // split the calculation to just 1 of 3 possible mutually exclusive ranges
  
   if (f < p->fInsJoin) // Inspiral range
   {
           double AmpIns = AmpPreFac * AmpInsAnsatz(f, powers_of_f, prefactors);
           return AmpIns;
   }
  
   if (f >= p->fMRDJoin) // MRD range
   {
           double AmpMRD = AmpPreFac * AmpMRDAnsatz(f, p);
           return AmpMRD;
   }
  
   //    Intermediate range
   double AmpInt = AmpPreFac * AmpIntAnsatz(f, p);
   return AmpInt;
}

/**
 * Take the AmpInsAnsatz expression and compute the first derivative
 * with respect to frequency to get the expression below.
 */
 double DAmpInsAnsatz(double Mf, useful_powers_s *powers_of_Mf, IMRPhenomDAmplitudeCoefficients* p) {
   double eta = p->eta;
   double chi1 = p->chi1;
   double chi2 = p->chi2;
   double rho1 = p->inspiral[0];
   double rho2 = p->inspiral[1];
   double rho3 = p->inspiral[2];
  
   double chi12 = p->chi12;
   double chi22 = p->chi22;
   double eta2 = p->eta2;
   double eta3 = p->eta3;
   double Pi = LAL_PI;
   double Pi2 = powers_of_pi.two;
   double Seta = p->Seta;
   double SetaPlus1 = p->SetaPlus1;
  
    return ((-969 + 1804*eta)*powers_of_pi.two_thirds)/(1008.*powers_of_Mf->third)
    + ((chi1*(81*SetaPlus1 - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Pi)/48.
    + ((-27312085 - 10287648*chi22 - 10287648*chi12*SetaPlus1
    + 10287648*chi22*Seta + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta
    + 35371056*eta2)*powers_of_Mf->third*powers_of_pi.four_thirds)/6.096384e6
    + (5*powers_of_Mf->two_thirds*powers_of_pi.five_thirds*(chi2*(-285197*(-1 + Seta)
    + 4*(-91902 + 1579*Seta)*eta - 35632*eta2) + chi1*(285197*SetaPlus1
    - 4*(91902 + 1579*Seta)*eta - 35632*eta2) + 42840*(-1 + 4*eta)*Pi))/96768.
    - (Mf*Pi2*(-336*(-3248849057.0 + 2943675504*chi12 - 3339284256*chi1*chi2 + 2943675504*chi22)*eta2 - 324322727232*eta3
    - 7*(-177520268561 + 107414046432*chi22 + 107414046432*chi12*SetaPlus1 - 107414046432*chi22*Seta
    + 11087290368*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi)
    + 12*eta*(-545384828789.0 - 176491177632*chi1*chi2 + 202603761360*chi22 + 77616*chi12*(2610335 + 995766*Seta)
    - 77287373856*chi22*Seta + 5841690624*(chi1 + chi2)*Pi + 21384760320*Pi2)))/3.0042980352e10
    +(7.0/3.0)*powers_of_Mf->four_thirds*rho1 + (8.0/3.0)*powers_of_Mf->five_thirds*rho2 + 3.*powers_of_Mf->two*rho3;
}

/**
 * first frequency derivative of AmpMRDAnsatz
 */
 double DAmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
   double fRD = p->fRD;
   double fDM = p->fDM;
   double gamma1 = p->merger_ringdown[0];
   double gamma2 = p->merger_ringdown[1];
   double gamma3 = p->merger_ringdown[2];
  
   double fDMgamma3 = fDM * gamma3;
   double pow2_fDMgamma3 = Square(fDMgamma3);
   double fminfRD = f - fRD;
   double expfactor = exp(((fminfRD)*gamma2)/(fDMgamma3));
   double pow2pluspow2 = Square(fminfRD) + pow2_fDMgamma3;
  
   return ((-2*fDM*(fminfRD)*gamma3*gamma1) / pow2pluspow2 -
     (gamma2*gamma1)) / ( expfactor * (pow2pluspow2)) ;
}


/**
 * Formula to predict the final spin. Equation 3.6 arXiv:1508.07250
 * s defined around Equation 3.6.
 */
 double ZombFinalSpin0815_s(double eta, double s) {
   double eta2 = eta*eta;
   double eta3 = eta2*eta;
   double s2 = s*s;
   double s3 = s2*s;
  
/* FIXME: there are quite a few int's withouth a . in this file */
//FP: eta2, eta3 can be avoided
return eta*(3.4641016151377544 - 4.399247300629289*eta +
       9.397292189321194*eta2 - 13.180949901606242*eta3 +
       s*((1.0/eta - 0.0850917821418767 - 5.837029316602263*eta) +
       (0.1014665242971878 - 2.0967746996832157*eta)*s +
       (-1.3546806617824356 + 4.108962025369336*eta)*s2 +
       (-0.8676969352555539 + 2.064046835273906*eta)*s3));
}

 double ZombFinalSpin0815(double eta, double chi1, double chi2) {
   // Convention m1 >= m2
   double Seta = sqrt(1.0 - 4.0*eta);
   double m1 = 0.5 * (1.0 + Seta);
   double m2 = 0.5 * (1.0 - Seta);
   double m1s = m1*m1;
   double m2s = m2*m2;
   // s defined around Equation 3.6 arXiv:1508.07250
   double s = (m1s * chi1 + m2s * chi2);
   return ZombFinalSpin0815_s(eta, s);
}


 int init_useful_powers(useful_powers_s *p, double number)
{
   XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
   XLAL_CHECK(number >= 0, XLAL_EDOM, "number must be non-negative");
  
   // consider changing pow(x,1/6.0) to cbrt(x) and sqrt(x) - might be faster
   double sixth = pow(number, 1.0 / 6.0);
   p->third = sixth * sixth;
   //p->third = cbrt(number);
   p->two_thirds = p->third * p->third;
   p->four_thirds = number * p->third;
   p->five_thirds = p->four_thirds * p->third;
   p->two = number * number;
   p->seven_thirds = p->third * p->two;
   p->eight_thirds = p->two_thirds * p->two;
   p->inv = 1. / number;
   double m_sixth = 1.0 / sixth;
   p->m_seven_sixths = p->inv * m_sixth;
   p->m_third = m_sixth * m_sixth;
   p->m_two_thirds = p->m_third * p->m_third;
   p->m_five_thirds = p->inv * p->m_two_thirds;
  
   return XLAL_SUCCESS;
}

 double fmaxCalc(IMRPhenomDAmplitudeCoefficients* p) {
   double fRD = p->fRD;
   double fDM = p->fDM;
   double gamma2 = p->merger_ringdown[1];
   double gamma3 = p->merger_ringdown[2];
  
   // NOTE: There's a problem with this expression from the paper becoming imaginary if gamma2>=1
   // Fix: if gamma2 >= 1 then set the Square root term to zero.
   if (!(gamma2 > 1))
     return fabs(fRD + (fDM*(-1 + sqrt(1 - Square(gamma2)))*gamma3)/gamma2);
   else
     return fabs(fRD + (-fDM*gamma3)/gamma2);
}

/******************************* Amplitude functions *******************************/
 
/**
 * amplitude scaling factor defined by eq. 17 in 1508.07253
 */
 double amp0Func(double eta) {
   return sqrt(2.0/3.0*eta)*PI_M_SIXTH;
}

 int init_amp_ins_prefactors(AmpInsPrefactors * prefactors, IMRPhenomDAmplitudeCoefficients* p)
{
         XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
         XLAL_CHECK(0 != prefactors, XLAL_EFAULT, "prefactors is NULL");
  
         double eta = p->eta;
  
         prefactors->amp0 = amp0Func(eta);
  
         double chi1 = p->chi1;
         double chi2 = p->chi2;
  
         double chi12 = p->chi12;
         double chi22 = p->chi22;
         double eta2 = p->eta2;
         double eta3 = p->eta3;
  
         double Pi = LAL_PI;
         double Pi2 = powers_of_pi.two;
         double Seta = p->Seta;
   double SetaPlus1 = p->SetaPlus1;
  
         prefactors->two_thirds = ((-969 + 1804*eta)*powers_of_pi.two_thirds)/672.;
         prefactors->one = ((chi1*(81*SetaPlus1 - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Pi)/48.;
         prefactors->four_thirds = (     (-27312085.0 - 10287648*chi22 - 10287648*chi12*SetaPlus1 + 10287648*chi22*Seta
                                                                  + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta
                                                                  + 35371056*eta2
                                                                  )
                                                         * powers_of_pi.four_thirds) / 8.128512e6;
         prefactors->five_thirds = (powers_of_pi.five_thirds * (chi2*(-285197*(-1 + Seta) + 4*(-91902 + 1579*Seta)*eta - 35632*eta2)
                                                                                                                         + chi1*(285197*SetaPlus1 - 4*(91902 + 1579*Seta)*eta - 35632*eta2)
                                                                                                                         + 42840*(-1.0 + 4*eta)*Pi
                                                                                                                         )
                                                                 ) / 32256.;
         prefactors->two = - (Pi2*(-336*(-3248849057.0 + 2943675504*chi12 - 3339284256*chi1*chi2 + 2943675504*chi22)*eta2
                                                           - 324322727232*eta3
                                                           - 7*(-177520268561 + 107414046432*chi22 + 107414046432*chi12*SetaPlus1
                                                                         - 107414046432*chi22*Seta + 11087290368*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi
                                                                         )
                                                           + 12*eta*(-545384828789 - 176491177632*chi1*chi2 + 202603761360*chi22
                                                                                 + 77616*chi12*(2610335 + 995766*Seta) - 77287373856*chi22*Seta
                                                                                 + 5841690624*(chi1 + chi2)*Pi + 21384760320*Pi2
                                                                                 )
                                                                 )
                                                 )/6.0085960704e10;
         prefactors->seven_thirds = p->inspiral[0];
         prefactors->eight_thirds = p->inspiral[1];
         prefactors->three        = p->inspiral[2];
  
         return XLAL_SUCCESS;
}

void ComputeDeltasFromCollocation(IMRPhenomDAmplitudeCoefficients* p) {
    // Three evenly spaced collocation points in the interval [f1,f3].
    double f1 = AMP_fJoin_INS;
    double f3 = p->fmaxCalc;
    double dfx = 0.5*(f3 - f1);
    double f2 = f1 + dfx;

    useful_powers_s powers_of_f1;
    int status = init_useful_powers(&powers_of_f1, f1);
    XLAL_CHECK_VOID ( status == XLAL_SUCCESS, XLAL_EFUNC, "Failed to initialize useful powers of f1.");

    AmpInsPrefactors prefactors;
    status = init_amp_ins_prefactors(&prefactors, p);
    XLAL_CHECK_VOID ( status == XLAL_SUCCESS, XLAL_EFUNC, "Failed to initialize amplitude prefactors for inspiral range.");

    // v1 is inspiral model evaluated at f1
    // d1 is derivative of inspiral model evaluated at f1
    double v1 = AmpInsAnsatz(f1, &powers_of_f1, &prefactors);
    double d1 = DAmpInsAnsatz(f1, &powers_of_f1, p);

    // v3 is merger-ringdown model evaluated at f3
    // d2 is derivative of merger-ringdown model evaluated at f3
    double v3 = AmpMRDAnsatz(f3, p);
    double d2 = DAmpMRDAnsatz(f3, p);

    // v2 is the value of the amplitude evaluated at f2
    // they come from the fit of the collocation points in the intermediate region

    double v2 = 
       calcCoefficient(
           AMPLITUDE_INTERMEDIATE_COLLOCATION_FIT_COEFFICIENT_TERMS, 
           p->eta, 
           p->eta2, 
           (-1.0 + p->chi)
       );  

    p->f1 = f1;
    p->f2 = f2;
    p->f3 = f3;
    p->v1 = v1;
    p->v2 = v2;
    p->v3 = v3;
    p->d1 = d1;
    p->d2 = d2;

    // Now compute the delta_i's from the collocation coefficients
    // Precompute common quantities here and pass along to delta functions.
    const double f1_1 = f1;
    const double f1_2 = f1*f1;
    const double f1_3 = f1*f1_2;
    const double f1_4 = f1*f1_3;
    const double f1_5 = f1*f1_4;
    
    const double f2_1 = f2;
    const double f2_2 = f2*f2;
    const double f2_3 = f2*f2_2;
    const double f2_4 = f2*f2_3;
    
    const double f3_1 = f3;
    const double f3_2 = f3*f3;
    const double f3_3 = f3*f3_2;
    const double f3_4 = f3*f3_3;
    const double f3_5 = f3*f3_4;
    
    // Use inner term macros defined in phenomd_data.h to abstract maths:
    const double inner_term_results[] =
    {
        DELTA_0_TERMS,
        DELTA_1_TERMS,
        DELTA_2_TERMS,
        DELTA_3_TERMS,
        DELTA_4_TERMS
    };
    
    const double divisor = (Square(p->f1 - p->f2)*Cube(p->f1 - p->f3)*Square(p->f3 - p->f2));
    for (int32_t index = 0; index < NUM_DELTA_TERMS; index++)
    {
        p->intermediate[index] = -1.0*inner_term_results[index] / divisor;
    }     
}

void ZombComputeIMRPhenomDAmplitudeCoefficients(IMRPhenomDAmplitudeCoefficients *p, double eta, double chi1, double chi2, double finspin) {

    p->eta = eta;
    p->etaInv = 1./eta;
    p->chi1 = chi1;
    p->chi2 = chi2;
    p->chi12 = chi1*chi1;
    p->chi22 = chi2*chi2;
    double eta2 = eta*eta;
    p->eta2 = eta2;
    p->eta3 = eta*eta2;
    double Seta = sqrt(1.0 - 4.0*eta);
    p->Seta = Seta;
    p->SetaPlus1 = 1.0 + Seta;

    p->q = 0.5 * (1.0 + Seta - 2.0*eta) * p->etaInv;
    p->chi = chiPN(Seta, eta, chi1, chi2);
    double xi = -1.0 + p->chi;

    p->fRD = Zombfring(eta, chi1, chi2, finspin);
    p->fDM = Zombfdamp(eta, chi1, chi2, finspin);
    
    // Compute gamma_i's, rho_i's first then delta_i's:
    calcCoefficients(
        AMPLITUDE_MERGER_RINGDOWN_COEFFICIENT_TERMS,
        NUM_AMPLITUDE_MERGER_RINGDOWN_COEFFICIENTS,
        eta, eta2, xi,
        p->merger_ringdown
    );
    calcCoefficients(
        AMPLITUDE_INSPIRAL_COEFFICIENT_TERMS,
        NUM_AMPLITUDE_INSPIRAL_COEFFICIENTS,
        eta, eta2, xi,
        p->inspiral
    );

    p->fmaxCalc = fmaxCalc(p);

    // compute delta_i's
    ComputeDeltasFromCollocation(p);
}

/**
 * This function computes the IMR phase given phenom coefficients.
 * Defined in VIII. Full IMR Waveforms arXiv:1508.07253
 * Rholm was added when IMRPhenomHM (high mode) was added.
 * Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
 * to 1. and PhenomD is recovered.
 * Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
 * Again, when Taulm = 1.0 then PhenomD is recovered.
 */

/**
 * Ansatz for the inspiral phase.
 * We call the LAL TF2 coefficients here.
 * The exact values of the coefficients used are given
 * as comments in the top of this file
 * Defined by Equation 27 and 28 arXiv:1508.07253
 */
 double PhiInsAnsatzInt(double Mf, useful_powers_s *powers_of_Mf, PhiInsPrefactors *prefactors, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn)
{
   XLAL_CHECK(0 != pn, XLAL_EFAULT, "pn is NULL");
       
   // Assemble PN phasing series
   const double v = powers_of_Mf->third * powers_of_pi.third;
   const double logv = log(v);
  
   double phasing = prefactors->initial_phasing;
  
   phasing += prefactors->two_thirds * powers_of_Mf->two_thirds;
   phasing += prefactors->third * powers_of_Mf->third;
   phasing += prefactors->third_with_logv * logv * powers_of_Mf->third;
   phasing += prefactors->logv * logv;
   phasing += prefactors->minus_third * powers_of_Mf->m_third;
   phasing += prefactors->minus_two_thirds * powers_of_Mf->m_two_thirds;
   phasing += prefactors->minus_one * powers_of_Mf->inv;
   phasing += prefactors->minus_four_thirds / powers_of_Mf->four_thirds;
   phasing += prefactors->minus_five_thirds * powers_of_Mf->m_five_thirds; // * v^0
  
   // Now add higher order terms that were calibrated for PhenomD
   phasing += ( prefactors->one * Mf + prefactors->four_thirds * powers_of_Mf->four_thirds
                            + prefactors->five_thirds * powers_of_Mf->five_thirds
                            + prefactors->two * powers_of_Mf->two
         ) * p->etaInv;
     
   return phasing;
}

 double DPhiIntAnsatz(double Mf, IMRPhenomDPhaseCoefficients *p) {
   return (p->beta1 + p->beta3/Quart(Mf) + p->beta2/Mf) * p->etaInv;
}

/**
 * ansatz for the intermediate phase defined by Equation 16 arXiv:1508.07253
 */
 double PhiIntAnsatz(double Mf, IMRPhenomDPhaseCoefficients *p) {
   // 1./eta in paper omitted and put in when need in the functions:
   // ComputeIMRPhenDPhaseConnectionCoefficients
   // IMRPhenDPhase
   return  p->beta1*Mf - p->beta3/(3.*Cube(Mf)) + p->beta2*log(Mf);
}

/**
 * Ansatz for the merger-ringdown phase Equation 14 arXiv:1508.07253
 * Rholm was added when IMRPhenomHM (high mode) was added.
 * Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
 * to 1. and PhenomD is recovered.
 * Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
 * Again, when Taulm = 1.0 then PhenomD is recovered.
 */
double PhiMRDAnsatzInt(double f, IMRPhenomDPhaseCoefficients *p, double Rholm, double Taulm)
{
   double sqrootf = sqrt(f);
   double fpow1_5 = f * sqrootf;
   // check if this is any faster: 2 sqrts instead of one pow(x,0.75)
   double fpow0_75 = sqrt(fpow1_5); // pow(f,0.75);
  
   return -(p->alpha2/f)
                  + (4.0/3.0) * (p->alpha3 * fpow0_75)
                  + p->alpha1 * f
                  + p->alpha4 * Rholm * atan((f - p->alpha5 * p->fRD) / (Rholm * p->fDM * Taulm));
}

 double IMRPhenDPhase(double f, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn, useful_powers_s *powers_of_f, PhiInsPrefactors *prefactors, double Rholm, double Taulm)
{
   // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
   // The inspiral, intermendiate and merger-ringdown phase parts
  
   // split the calculation to just 1 of 3 possible mutually exclusive ranges
   if (f < p->fInsJoin)        // Inspiral range
   {
           double PhiIns = PhiInsAnsatzInt(f, powers_of_f, prefactors, p, pn);
       
           return PhiIns;
   }
  
   if (f >= p->fMRDJoin) // MRD range
   {
           double PhiMRD = p->etaInv * PhiMRDAnsatzInt(f, p, Rholm, Taulm) + p->C1MRD + p->C2MRD * f;
       

           return PhiMRD;
   }
  
   //    Intermediate range
   double PhiInt = p->etaInv * PhiIntAnsatz(f, p) + p->C1Int + p->C2Int * f;
   return PhiInt;
}

 void ComputeIMRPhenomDPhaseCoefficients(IMRPhenomDPhaseCoefficients *p, double eta, double chi1, double chi2, double finspin, LALDict *extraParams) {
  
   // Convention m1 >= m2
   p->eta = eta;
   p->etaInv = 1./eta;
   p->chi1 = chi1;
   p->chi2 = chi2;
   double eta2 = eta*eta;
   p->eta2 = eta2;
   p->Seta = sqrt(1.0 - 4.0*eta);
  
   p->q = 0.5*(1.0 + p->Seta - 2.0*eta)*p->etaInv;
   p->chi = chiPN(p->Seta, eta, chi1, chi2);
   double xi = -1.0 + p->chi;
  
   p->sigma1 = calcCoefficient(PHASE_INSPIRAL_COEFFICIENT_TERMS[0], eta, eta2, xi);
   p->sigma2 = calcCoefficient(PHASE_INSPIRAL_COEFFICIENT_TERMS[1], eta, eta2, xi);
   p->sigma3 = calcCoefficient(PHASE_INSPIRAL_COEFFICIENT_TERMS[2], eta, eta2, xi);
   p->sigma4 = calcCoefficient(PHASE_INSPIRAL_COEFFICIENT_TERMS[3], eta, eta2, xi);
  
   p->beta1 = calcCoefficient(PHASE_INTERMEDIATE_COEFFICIENT_TERMS[0], eta, eta2, xi);
   p->beta2 = calcCoefficient(PHASE_INTERMEDIATE_COEFFICIENT_TERMS[1], eta, eta2, xi);
   p->beta3 = calcCoefficient(PHASE_INTERMEDIATE_COEFFICIENT_TERMS[2], eta, eta2, xi);
  
   p->alpha1 = calcCoefficient(PHASE_RINGDOWN_COEFFICIENT_TERMS[0], eta, eta2, xi);
   p->alpha2 = calcCoefficient(PHASE_RINGDOWN_COEFFICIENT_TERMS[1], eta, eta2, xi);
   p->alpha3 = calcCoefficient(PHASE_RINGDOWN_COEFFICIENT_TERMS[2], eta, eta2, xi);
   p->alpha4 = calcCoefficient(PHASE_RINGDOWN_COEFFICIENT_TERMS[3], eta, eta2, xi);
   p->alpha5 = calcCoefficient(PHASE_RINGDOWN_COEFFICIENT_TERMS[4], eta, eta2, xi);
  
   p->fRD = Zombfring(eta, chi1, chi2, finspin);
   p->fDM = Zombfdamp(eta, chi1, chi2, finspin);
  
   p->sigma1*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDSigma1(extraParams));
   p->sigma2*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDSigma2(extraParams));
   p->sigma3*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDSigma3(extraParams));
   p->sigma4*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDSigma4(extraParams));
   p->beta1*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDBeta1(extraParams));
   p->beta2*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDBeta2(extraParams));
   p->beta3*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDBeta3(extraParams));
   p->alpha1*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha1(extraParams));
   p->alpha2*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha2(extraParams));
   p->alpha3*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha3(extraParams));
   p->alpha4*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha4(extraParams));
   p->alpha5*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDAlpha5(extraParams));
  
}

 int init_phi_ins_prefactors(PhiInsPrefactors * prefactors, IMRPhenomDPhaseCoefficients* p, PNPhasingSeries *pn)
{
         XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
         XLAL_CHECK(0 != prefactors, XLAL_EFAULT, "prefactors is NULL");
  
         double sigma1 = p->sigma1;
         double sigma2 = p->sigma2;
         double sigma3 = p->sigma3;
         double sigma4 = p->sigma4;
  
   // PN phasing series
         prefactors->initial_phasing = pn->v[5] - LAL_PI_4;
         prefactors->two_thirds = pn->v[7] * powers_of_pi.two_thirds;
         prefactors->third = pn->v[6] * powers_of_pi.third;
         prefactors->third_with_logv = pn->vlogv[6] * powers_of_pi.third;
         prefactors->logv = pn->vlogv[5];
         prefactors->minus_third = pn->v[4] * powers_of_pi.m_third;
         prefactors->minus_two_thirds = pn->v[3] * powers_of_pi.m_two_thirds;
         prefactors->minus_one = pn->v[2] * powers_of_pi.inv;
   prefactors->minus_four_thirds = pn->v[1] / powers_of_pi.four_thirds;
   prefactors->minus_five_thirds = pn->v[0] * powers_of_pi.m_five_thirds; // * v^0
  
   // higher order terms that were calibrated for PhenomD
         prefactors->one = sigma1;
         prefactors->four_thirds = sigma2 * 0.75;
         prefactors->five_thirds = sigma3 * 0.6;
         prefactors->two = sigma4 * 0.5;
  
         return XLAL_SUCCESS;
}

/**
 * temporary instance of DPhiIntAnsatz used when computing
 * coefficients to make the phase C(1) continuous between regions.
 */
 double DPhiIntTemp(double ff, IMRPhenomDPhaseCoefficients *p) {
   double etaInv = p->etaInv;
   double beta1 = p->beta1;
   double beta2 = p->beta2;
   double beta3 = p->beta3;
   double C2Int = p->C2Int;
  
   return C2Int + (beta1 + beta3/Quart(ff) + beta2/ff)*etaInv;
}

/**
 * First frequency derivative of PhiInsAnsatzInt
 */
 double DPhiInsAnsatzInt(double Mf, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn) {
   double sigma1 = p->sigma1;
   double sigma2 = p->sigma2;
   double sigma3 = p->sigma3;
   double sigma4 = p->sigma4;
   double Pi = LAL_PI;
  
   // Assemble PN phasing series
   const double v = cbrt(Pi*Mf);
   const double logv = log(v);
   const double v2 = v * v;
   const double v3 = v * v2;
   const double v4 = v * v3;
   const double v5 = v * v4;
   const double v6 = v * v5;
   const double v7 = v * v6;
   const double v8 = v * v7;
  
   // Apply the correct prefactors to LAL phase coefficients to get the
   // phase derivative dphi / dMf = dphi/dv * dv/dMf
   double Dphasing = 0.0;
   Dphasing += 2.0 * pn->v[7] * v7;
   Dphasing += (pn->v[6] + pn->vlogv[6] * (1.0 + logv)) * v6;
   Dphasing += pn->vlogv[5] * v5;
   Dphasing += -1.0 * pn->v[4] * v4;
   Dphasing += -2.0 * pn->v[3] * v3;
   Dphasing += -3.0 * pn->v[2] * v2;
   Dphasing += -4.0 * pn->v[1] * v;
   Dphasing += -5.0 * pn->v[0];
   Dphasing /= v8 * 3.0;
   Dphasing *= Pi;
  
   // Now add higher order terms that were calibrated for PhenomD
   Dphasing += (
           sigma1
         + sigma2 * v * powers_of_pi.m_third
         + sigma3 * v2 * powers_of_pi.m_two_thirds
         + (sigma4*powers_of_pi.inv) * v3
         ) * p->etaInv;
  
   return Dphasing;
}

/**
 * First frequency derivative of PhiMRDAnsatzInt
 * Rholm was added when IMRPhenomHM (high mode) was added.
 * Rholm = fRD22/fRDlm. For PhenomD (only (l,m)=(2,2)) this is just equal
 * to 1. and PhenomD is recovered.
 * Taulm = fDMlm/fDM22. Ratio of ringdown damping times.
 * Again, when Taulm = 1.0 then PhenomD is recovered.
 */
 double DPhiMRD(double f, IMRPhenomDPhaseCoefficients *p, double Rholm, double Taulm) {

   return ( p->alpha1 + p->alpha2/Square(f) + p->alpha3/pow(f,0.25)+ p->alpha4/(p->fDM * Taulm * (1 + Square(f - p->alpha5 * p->fRD)/(Square(p->fDM * Taulm * Rholm)))) ) * p->etaInv;
}

 void ComputeIMRPhenDPhaseConnectionCoefficients(IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn, PhiInsPrefactors *prefactors, double Rholm, double Taulm)
{
   double etaInv = p->etaInv;
  
   // Transition frequencies
   // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
   p->fInsJoin=PHI_fJoin_INS;
   p->fMRDJoin=0.5*p->fRD;
  
   // Compute C1Int and C2Int coeffs
   // Equations to solve for to get C(1) continuous join
   // PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
   // Joining at fInsJoin
   // PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
   // PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
   double DPhiIns = DPhiInsAnsatzInt(PHI_fJoin_INS, p, pn);
   double DPhiInt = DPhiIntAnsatz(PHI_fJoin_INS, p);
   p->C2Int = DPhiIns - DPhiInt;
  
   useful_powers_s powers_of_fInsJoin;
   init_useful_powers(&powers_of_fInsJoin, PHI_fJoin_INS);
   p->C1Int = PhiInsAnsatzInt(PHI_fJoin_INS, &powers_of_fInsJoin, prefactors, p, pn)
     - etaInv * PhiIntAnsatz(PHI_fJoin_INS, p) - p->C2Int * PHI_fJoin_INS;
  
   // Compute C1MRD and C2MRD coeffs
   // Equations to solve for to get C(1) continuous join
   // PhiInsInt (f)  =   PhiMRD (f) + C1MRD + C2MRD f
   // Joining at fMRDJoin
   // Where \[Phi]InsInt(f) is the \[Phi]Ins+\[Phi]Int joined function
   // PhiInsInt (fMRDJoin)  =   PhiMRD (fMRDJoin) + C1MRD + C2MRD fMRDJoin
   // PhiInsInt'(fMRDJoin)  =   PhiMRD'(fMRDJoin) + C2MRD
   // temporary Intermediate Phase function to Join up the Merger-Ringdown
   double PhiIntTempVal = etaInv * PhiIntAnsatz(p->fMRDJoin, p) + p->C1Int + p->C2Int*p->fMRDJoin;
   double DPhiIntTempVal = DPhiIntTemp(p->fMRDJoin, p);
   double DPhiMRDVal = DPhiMRD(p->fMRDJoin, p, Rholm, Taulm);
        
   p->C2MRD = DPhiIntTempVal - DPhiMRDVal;
   p->C1MRD = PhiIntTempVal - etaInv * PhiMRDAnsatzInt(p->fMRDJoin, p, Rholm, Taulm) - p->C2MRD*p->fMRDJoin;
}