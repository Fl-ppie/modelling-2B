# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:19:54 2024

@author: Fardeen
To do:
    
"""

import Galaxies_and_Dynamics as GD
import numpy as np

# Constants
M_sun = 1.989e30  # Solar mass in kg
kpc_to_m = 3.086e19  # Kiloparsec to meters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
a0 = 1.2e-10  # MOND acceleration threshold (m/s^2)

"""
Running the simulation
"""
# Simulation parameters
dt = 1e13  # Time step (seconds) - smaller to improve stability
total_time = 1e19  # Total time (in seconds) ~58 days
steps = int(total_time / dt)

# Run simulation
def run_simulation(galaxies, mode="simple"):
    positions = {i: [] for i in range(len(galaxies))}
    for step in range(steps):
        if step % 100 == 0:
            print(f"Step {step}/{steps}")
        GD.runge_kutta4(galaxies, dt, mode=mode)
        for i, galaxy in enumerate(galaxies):
            positions[i].append(galaxy.position.copy())
    for i in range(len(galaxies)):
        positions[i] = np.array(positions[i])
    return positions