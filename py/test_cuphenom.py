from cuphenom import generatePhenom
import numpy as np

from bokeh.layouts import gridplot
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource

from tqdm import tqdm

import time
import random

def test_generatePhenom():
    # Define reasonable input parameters for an average gravitational wave
    approximant_enum = 1
    mass_1_msun = 30
    mass_2_msun = 30
    sample_rate_hertz = 4096
    duration_seconds = 4.0
    inclination_radians = 1.0
    distance_mpc = 100.0
    reference_orbital_phase_in = 0.0
    ascending_node_longitude = 100.0
    eccentricity = 0.0
    mean_periastron_anomaly = 0.0
    spin_1_in = [0.0, 0.0, 0.0]
    spin_2_in = [0.0, 0.0, 0.0]

    # Call generatePhenom function
    result = generatePhenom(
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
    )  
    result = np.array(result.tolist())
    
    # Assume that result corresponds to pairs of gravitational wave polarisations plus and cross
    times = np.arange(0, duration_seconds, 1/sample_rate_hertz)
    polarisations_plus = result[:, 0]
    polarisations_cross = result[:, 1]

    # Prepare data for bokeh
    data = ColumnDataSource(data=dict(
        time=times,
        plus=polarisations_plus,
        cross=polarisations_cross
    ))

    # Create a new plot
    p = figure(title="Gravitational Wave Polarisations", x_axis_label='Time (s)', y_axis_label='Strain')

    # Add polarisation traces
    p.line('time', 'plus', source=data, legend_label="Plus Polarisation", line_color="blue")
    p.line('time', 'cross', source=data, legend_label="Cross Polarisation", line_color="red")

    # Move the legend to the upper left corner
    p.legend.location = "top_left"

    # Output to static HTML file
    output_file("gravitational_wave_polarisations.html")

    # Show the results
    show(p)

def speed_test_generatePhenom(num_tests=100):

    # Prepare data storage for bokeh
    times = []
    runtimes = []

    for i in tqdm(range(num_tests)):
        # Define random input parameters
        approximant_enum = 0 # random choice between 1, 2, 3
        mass_1_msun = random.uniform(10, 50)  # random float between 1 and 30
        mass_2_msun = random.uniform(10, 50)
        sample_rate_hertz = 4096  # Keep these parameters fixed for now
        duration_seconds = 4.0
        inclination_radians = random.uniform(0, np.pi)
        distance_mpc = random.uniform(10.0, 1000.0)
        reference_orbital_phase_in = random.uniform(0, np.pi*2)
        ascending_node_longitude = random.uniform(0, np.pi*2)
        eccentricity = random.uniform(0, 0.1)
        mean_periastron_anomaly = random.uniform(0, np.pi*2)
        spin_1_in = [random.uniform(-0.5, 0.5) for _ in range(3)]
        spin_2_in = [random.uniform(-0.5, 0.5) for _ in range(3)]

        # Start the timer
        start_time = time.time()

        # Call generatePhenom function
        result = generatePhenom(
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
        )

        # Record the runtime
        runtime = time.time() - start_time
        times.append(start_time)
        runtimes.append(runtime)
        
    print("Runtimes", np.sum(runtimes))

    # Prepare data for bokeh
    data = ColumnDataSource(data=dict(
        time=times,
        runtime=runtimes,
    ))

    # Create a new plot
    p = figure(title="GeneratePhenom Run Time", x_axis_label='Start Time (s)', y_axis_label='Run Time (s)')

    # Add runtime traces
    p.circle('time', 'runtime', source=data, legend_label="Runtime", line_color="blue")

    # Move the legend to the upper left corner
    p.legend.location = "top_left"

    # Output to static HTML file
    output_file("generatePhenom_runtimes.html")

    # Show the results
    show(p)

# Call the test function
test_generatePhenom()
#speed_test_generatePhenom(num_tests=1000000)