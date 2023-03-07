#ifndef PHENOM_D_FUNCTIONS_H
#define PHENOM_D_FUNCTIONS_H

complex_waveform_axes_s generatePhenomD(
    const temporal_properties_s *temporal_properties,
    const int32_t                num_strain_axis_samples
    );

int32_t sumPhenomDFrequencies(
          complex_waveform_axes_s         waveform_axes_fd,
    const float                           inclination,
    const float                           total_mass_seconds,
    const amplitude_coefficients_s        amplitude_coefficients,
    const amplitude_inspiral_prefactors_s amplitude_prefactors,
    const phase_coefficients_s            phase_coefficients, 
    const phase_inspiral_prefactors_s     phase_prefactors, 
    const int32_t                         offset,
    const float                           phase_shift,
    const float                           amp0,
    const float                           reference_mass_frequency,
    const float                           phi_precalc
    );

m_complex_waveform_axes_s initPhenomDWaveformAxes(
          temporal_properties_s *temporal_properties,
          system_properties_s   *system_properties,
    const int32_t                num_waveforms
    );

#endif