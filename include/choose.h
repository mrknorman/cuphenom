/** Helper routines for XLALSimInspiralTD(): performs conditioning of a FD waveform and transforms it to TD */


int InspiralTDFromFD_old(
     REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
     REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
     REAL8 m1,                                   /**< mass of companion 1 (kg) */
     REAL8 m2,                                   /**< mass of companion 2 (kg) */
     REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
     REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
     REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
     REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
     REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
     REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
     REAL8 distance,                             /**< distance of source (m) */
     REAL8 inclination,                          /**< inclination of source (rad) */
     REAL8 phiRef,                               /**< reference orbital phase (rad) */
     REAL8 longAscNodes,                         /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
     REAL8 eccentricity,                         /**< eccentrocity at reference epoch */
     REAL8 meanPerAno,                           /**< mean anomaly of periastron */
     REAL8 deltaT,                               /**< sampling interval (s) */
     REAL8 f_min,                                /**< starting GW frequency (Hz) */
     REAL8 f_ref,                                /**< reference GW frequency (Hz) */
     LALDict *LALparams,                         /**< LAL dictionary containing accessory parameters */
     Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
)
{
     COMPLEX16FrequencySeries *hptilde = NULL;
     COMPLEX16FrequencySeries *hctilde = NULL;
     void *plan;
     double tshift;
     const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
     const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
     double original_f_min = f_min; /* f_min might be overwritten below, so keep original value */
     double f_max = 0.5 / deltaT;
     double tchirp, tmerge, textra;
     double fstart;
     double s;
     int retval;
    
     /* upper bound on the chirp time starting at f_min */
     tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, S1z, S2z);
  
     /* upper bound on the final black hole spin */
     s = XLALSimInspiralFinalBlackHoleSpinBound(S1z, S2z);
  
     /* upper bound on the final plunge, merger, and ringdown time */
     tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);
  
     /* extra time to include for all waveforms to take care of situations
     * where the frequency is close to merger (and is sweeping rapidly):
     * this is a few cycles at the low frequency */
     textra = extra_cycles / f_min;
  
     /* generate the conditioned waveform in the frequency domain */
     /* note: redshift factor has already been applied above */
     /* set deltaF = 0 to get a small enough resolution */
     retval = XLALSimInspiralFD(&hptilde, &hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, 0.0, f_min, f_max, f_ref, LALparams, approximant);
     if (retval < 0)
         XLAL_ERROR(XLAL_EFUNC);
  
     /* we want to make sure that this waveform will give something
     * sensible if it is later transformed into the time domain:
     * to avoid the end of the waveform wrapping around to the beginning,
     * we shift waveform backwards in time and compensate for this
     * shift by adjusting the epoch -- note that XLALSimInspiralFD
     * guarantees that there is extra padding to do this */
     tshift = round(textra / deltaT) * deltaT; /* integer number of samples */
     for (size_t k = 0; k < hptilde->data->length; ++k) {
         double complex phasefac = cexp(2.0 * M_PI * I * k * hptilde->deltaF * tshift);
         hptilde->data->data[k] *= phasefac;
         hctilde->data->data[k] *= phasefac;
     }
     XLALGPSAdd(&hptilde->epoch, tshift);
     XLALGPSAdd(&hctilde->epoch, tshift);
	
	const LALUnit lalStrainUnit = { 0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
  
     /* transform the waveform into the time domain */
     size_t chirplen = 2 * (hptilde->data->length - 1);
     *hplus = XLALCreateREAL8TimeSeries("H_PLUS", &hptilde->epoch, 0.0, deltaT, &lalStrainUnit, chirplen);
     *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &hctilde->epoch, 0.0, deltaT, &lalStrainUnit, chirplen);
     plan = XLALCreateReverseREAL8FFTPlan(chirplen, 0);
     if (!(*hplus) || !(*hcross) || !plan) {
         XLALDestroyCOMPLEX16FrequencySeries(hptilde);
         XLALDestroyCOMPLEX16FrequencySeries(hctilde);
         XLALDestroyREAL8TimeSeries(*hcross);
         XLALDestroyREAL8TimeSeries(*hplus);
         XLALDestroyREAL8FFTPlan(plan);
         XLAL_ERROR(XLAL_EFUNC);
     }
     XLALREAL8FreqTimeFFT(*hplus, hptilde, plan);
     XLALREAL8FreqTimeFFT(*hcross, hctilde, plan);
  
     /* apply time domain filter at original f_min */
     XLALHighPassREAL8TimeSeries(*hplus, original_f_min, 0.99, 8);
     XLALHighPassREAL8TimeSeries(*hcross, original_f_min, 0.99, 8);
  
     /* compute how long a chirp we should have */
     /* revised estimate of chirp length from new start frequency */
     fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp, m1, m2);
     tchirp = XLALSimInspiralChirpTimeBound(fstart, m1, m2, S1z, S2z);
  
     /* total expected chirp length includes merger */
    chirplen = (size_t) round((tchirp + tmerge) / deltaT);
  
     /* amount to snip off at the end is tshift */
     size_t end =  (*hplus)->data->length - (size_t) round(tshift / deltaT);
  
     /* snip off extra time at beginning and at the end */
     XLALResizeREAL8TimeSeries(*hplus, end - chirplen, chirplen);
     XLALResizeREAL8TimeSeries(*hcross, end - chirplen, chirplen);
  
     /* clean up */
     XLALDestroyREAL8FFTPlan(plan);
     XLALDestroyCOMPLEX16FrequencySeries(hptilde);
     XLALDestroyCOMPLEX16FrequencySeries(hctilde);
  
     return 0;
}


static int XLALSimInspiralTDFromFD(
     REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
     REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
     REAL8 m1,                                   /**< mass of companion 1 (kg) */
     REAL8 m2,                                   /**< mass of companion 2 (kg) */
     REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
     REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
     REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
     REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
     REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
     REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
     REAL8 distance,                             /**< distance of source (m) */
     REAL8 inclination,                          /**< inclination of source (rad) */
     REAL8 phiRef,                               /**< reference orbital phase (rad) */
     REAL8 longAscNodes,                         /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
     REAL8 eccentricity,                         /**< eccentrocity at reference epoch */
     REAL8 meanPerAno,                           /**< mean anomaly of periastron */
     REAL8 deltaT,                               /**< sampling interval (s) */
     REAL8 f_min,                                /**< starting GW frequency (Hz) */
     REAL8 f_ref,                                /**< reference GW frequency (Hz) */
     LALDict *LALparams,                         /**< LAL dictionary containing accessory parameters */
     Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
)
{
     COMPLEX16FrequencySeries *hptilde = NULL;
     COMPLEX16FrequencySeries *hctilde = NULL;
     void *plan;
     size_t chirplen, end, k;
     double tshift;
     const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
     const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
     double original_f_min = f_min; /* f_min might be overwritten below, so keep original value */
     double f_max = 0.5 / deltaT;
     double tchirp, tmerge, textra;
     double fisco, fstart;
     double s;
     int retval;
  
     if (!XLALSimInspiralImplementedFDApproximants(approximant))
         XLAL_ERROR(XLAL_EINVAL, "Invalid approximant: not a FD approximant");
  
     /* apply redshift correction to dimensionful source-frame quantities */
     REAL8 z=XLALSimInspiralWaveformParamsLookupRedshift(LALparams);
     if (z != 0.0) {
         m1 *= (1.0 + z);
         m2 *= (1.0 + z);
         distance *= (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
     }
     /* set redshift to zero so we don't accidentally apply it again later */
     z=0.;
     if (LALparams)
       XLALSimInspiralWaveformParamsInsertRedshift(LALparams,z);
  
     /* if the requested low frequency is below the lowest Kerr ISCO
     * frequency then change it to that frequency */
     fisco = 1.0 / (pow(9.0, 1.5) * M_PI * (m1 + m2) * MASS_SUN_SECONDS / MASS_SUN_KILOGRAMS);
     if (f_min > fisco)
         f_min = fisco;
  
     /* upper bound on the chirp time starting at f_min */
     tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, S1z, S2z);
  
     /* upper bound on the final black hole spin */
     s = XLALSimInspiralFinalBlackHoleSpinBound(S1z, S2z);
  
     /* upper bound on the final plunge, merger, and ringdown time */
     tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);
  
     /* extra time to include for all waveforms to take care of situations
     * where the frequency is close to merger (and is sweeping rapidly):
     * this is a few cycles at the low frequency */
     textra = extra_cycles / f_min;
  
     /* generate the conditioned waveform in the frequency domain */
     /* note: redshift factor has already been applied above */
     /* set deltaF = 0 to get a small enough resolution */
     retval = XLALSimInspiralFD(&hptilde, &hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, 0.0, f_min, f_max, f_ref, LALparams, approximant);
     if (retval < 0)
         XLAL_ERROR(XLAL_EFUNC);
  
     /* we want to make sure that this waveform will give something
     * sensible if it is later transformed into the time domain:
     * to avoid the end of the waveform wrapping around to the beginning,
     * we shift waveform backwards in time and compensate for this
     * shift by adjusting the epoch -- note that XLALSimInspiralFD
     * guarantees that there is extra padding to do this */
     tshift = round(textra / deltaT) * deltaT; /* integer number of samples */
     for (k = 0; k < hptilde->data->length; ++k) {
         double complex phasefac = cexp(2.0 * M_PI * I * k * hptilde->deltaF * tshift);
         hptilde->data->data[k] *= phasefac;
         hctilde->data->data[k] *= phasefac;
     }
     XLALGPSAdd(&hptilde->epoch, tshift);
     XLALGPSAdd(&hctilde->epoch, tshift);
	  
     const LALUnit lalStrainUnit = { 0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };

     /* transform the waveform into the time domain */
     chirplen = 2 * (hptilde->data->length - 1);
     *hplus = XLALCreateREAL8TimeSeries("H_PLUS", &hptilde->epoch, 0.0, deltaT, &lalStrainUnit, chirplen);
     *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &hctilde->epoch, 0.0, deltaT, &lalStrainUnit, chirplen);
     plan = XLALCreateReverseREAL8FFTPlan(chirplen, 0);
     if (!(*hplus) || !(*hcross) || !plan) {
         XLALDestroyCOMPLEX16FrequencySeries(hptilde);
         XLALDestroyCOMPLEX16FrequencySeries(hctilde);
         XLALDestroyREAL8TimeSeries(*hcross);
         XLALDestroyREAL8TimeSeries(*hplus);
         XLALDestroyREAL8FFTPlan(plan);
         XLAL_ERROR(XLAL_EFUNC);
     }
     XLALREAL8FreqTimeFFT(*hplus, hptilde, plan);
     XLALREAL8FreqTimeFFT(*hcross, hctilde, plan);
  
     /* apply time domain filter at original f_min */
     XLALHighPassREAL8TimeSeries(*hplus, original_f_min, 0.99, 8);
     XLALHighPassREAL8TimeSeries(*hcross, original_f_min, 0.99, 8);
  
     /* compute how long a chirp we should have */
     /* revised estimate of chirp length from new start frequency */
     fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp, m1, m2);
     tchirp = XLALSimInspiralChirpTimeBound(fstart, m1, m2, S1z, S2z);
  
     /* total expected chirp length includes merger */
     chirplen = round((tchirp + tmerge) / deltaT);
  
     /* amount to snip off at the end is tshift */
     end = (*hplus)->data->length - round(tshift / deltaT);
  
     /* snip off extra time at beginning and at the end */
     XLALResizeREAL8TimeSeries(*hplus, end - chirplen, chirplen);
     XLALResizeREAL8TimeSeries(*hcross, end - chirplen, chirplen);
  
     /* clean up */
     XLALDestroyREAL8FFTPlan(plan);
     XLALDestroyCOMPLEX16FrequencySeries(hptilde);
     XLALDestroyCOMPLEX16FrequencySeries(hctilde);
  
     /* final tapering at the beginning and at the end to remove filter transients */
  
     /* waveform should terminate at a frequency >= Schwarzschild ISCO
     * so taper one cycle at this frequency at the end; should not make
     * any difference to IMR waveforms */
     fisco = 1.0 / (pow(6.0, 1.5) * M_PI * (m1 + m2) * MASS_SUN_SECONDS / MASS_SUN_KILOGRAMS);
     XLALSimInspiralTDConditionStage2(*hplus, *hcross, f_min, fisco);
  
     return 0;
}