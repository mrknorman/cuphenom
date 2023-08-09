import numpy as np
from ctypes import *
import os

# Load the library
current_dir = os.path.dirname(os.path.realpath(__file__))
lib = CDLL(f'{current_dir}/libphenom.so')

# Define the float2 structure
class Float2(Structure):
    _fields_ = [('x', c_float),
                ('y', c_float)]

lib.pythonWrapperPhenomD.argtypes = [
    c_int,  # approximant_enum
    c_float,  # mass_1_msun
    c_float,  # mass_2_msun
    c_float,  # sample_rate_hertz
    c_float,  # duration_seconds
    c_float,  # inclination_radians
    c_float,  # distance_mpc
    c_float,  # reference_orbital_phase_in
    c_float,  # ascending_node_longitude
    c_float,  # eccentricity
    c_float,  # mean_periastron_anomaly
    POINTER(c_float),  # spin_1_in
    POINTER(c_float)  # spin_2_in
]
lib.pythonWrapperPhenomD.restype = POINTER(c_float)

def generate_phenom(
	approximant_enum,
	mass_1_msun,
	mass_2_msun,
	sample_rate_hertz,
	duration_seconds,
	inclination_radians,
	distance_mpc,
	reference_orbital_phase_in,
	ascending_node_longitude,
	eccentricity,
	mean_periastron_anomaly,
	spin_1_in,
	spin_2_in
	):

    # Ensure input spins are float arrays of length 3
    assert len(spin_1_in) == len(spin_2_in) == 3

    # Convert numpy arrays to ctypes arrays
    spin_1_in_ctypes = (c_float * len(spin_1_in))(*spin_1_in)
    spin_2_in_ctypes = (c_float * len(spin_2_in))(*spin_2_in)

    # Call the C function
    result_ptr = lib.pythonWrapperPhenomD(int(approximant_enum),
                                          mass_1_msun,
                                          mass_2_msun,
                                          sample_rate_hertz,
                                          duration_seconds,
                                          inclination_radians,
                                          distance_mpc,
                                          reference_orbital_phase_in,
                                          ascending_node_longitude,
                                          eccentricity,
                                          mean_periastron_anomaly,
                                          spin_1_in_ctypes,
                                          spin_2_in_ctypes)

    # Since we don't know the size of the returned array, let's assume it's num_samples
    num_samples = int(np.floor(sample_rate_hertz * duration_seconds))
    
    result = np.ctypeslib.as_array(result_ptr, shape=(num_samples * 2,))
    result = np.column_stack((result[::2], result[1::2]))
        
    # Don't forget to free the memory if it's no longer needed
    # Your C library needs to provide a function for this; let's assume it's called "free"
    lib.free(result_ptr)

    return result
    
    