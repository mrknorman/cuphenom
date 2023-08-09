import numpy as np
from ctypes import *
import os

def to_ctypes(to_convert):
    return (c_float * len(to_convert))(*to_convert)

def ensure_list(item):
    return item if isinstance(item, (list, tuple, set, np.ndarray)) else [item]

# Load the library
current_dir = os.path.dirname(os.path.realpath(__file__))
lib = CDLL(f'{current_dir}/libphenom.so')

lib.pythonWrapperPhenomD.argtypes = \
[
    c_int,    # Num Waveforms
    c_float,  # sample_rate_hertz
    c_float,  # duration_seconds
    POINTER(c_float),  # mass_1_msun
    POINTER(c_float),  # mass_2_msun
    POINTER(c_float),  # inclination_radians
    POINTER(c_float),  # distance_mpc
    POINTER(c_float),  # reference_orbital_phase_in
    POINTER(c_float),  # ascending_node_longitude
    POINTER(c_float),  # eccentricity
    POINTER(c_float),  # mean_periastron_anomaly
    POINTER(c_float),  # spin_1_in
    POINTER(c_float)  # spin_2_in
]
lib.pythonWrapperPhenomD.restype = POINTER(c_float)

def generate_phenom_d(
    num_waveforms,
    sample_rate_hertz,
	duration_seconds,
	mass_1_msun,
	mass_2_msun,
	inclination_radians,
	distance_mpc,
	reference_orbital_phase_in,
	ascending_node_longitude,
	eccentricity,
	mean_periastron_anomaly,
	spin_1_in,
	spin_2_in
	):
    
    args = locals().copy()
    
    length_one_args = ['num_waveforms', 'sample_rate_hertz', 'duration_seconds']
    length_three_args = ['spin_1_in', 'spin_2_in']
    
    args = {name: ensure_list(val) for name, val in args.items()}
            
    # Check that all arguments have length num_waveforms
    args_to_check = {
        name: val for name, val in args.items()
        if name not in length_one_args + length_three_args
    }
    
    for arg_name, arg_value in args_to_check.items():
        assert len(arg_value) == num_waveforms, \
        f"{arg_name} does not have the expected length of {num_waveforms}, instead = {len(arg_value)}"
    
    # Ensure input spins are float arrays of length 3
    assert len(spin_1_in) == len(spin_2_in) == 3*num_waveforms
    
    # Convert all arguments except the first three to ctypes
    ctypes_args = {
        name: to_ctypes(val) for name, val in args.items()
        if name not in length_one_args
    }

    # Call the C function
    result_ptr = lib.pythonWrapperPhenomD(
        num_waveforms,
        sample_rate_hertz,
        duration_seconds,
        ctypes_args["mass_1_msun"],
        ctypes_args["mass_2_msun"],
        ctypes_args["inclination_radians"],
        ctypes_args["distance_mpc"],
        ctypes_args["reference_orbital_phase_in"],
        ctypes_args["ascending_node_longitude"],
        ctypes_args["eccentricity"],
        ctypes_args["mean_periastron_anomaly"],
        ctypes_args["spin_1_in"],
        ctypes_args["spin_2_in"]
    )

    num_samples = \
        int(np.floor(sample_rate_hertz * duration_seconds))*num_waveforms
    
    result = np.ctypeslib.as_array(result_ptr, shape=(num_samples * 2,))
    result = np.column_stack((result[::2], result[1::2]))
    
    lib.free(result_ptr)

    return result
    
    