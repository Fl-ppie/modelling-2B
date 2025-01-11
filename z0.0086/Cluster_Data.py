# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:33:51 2024

@author: Fardeen

To Do: 
    
"""

import Galaxies_and_Dynamics as GD
import numpy as np

# Constants
M_sun = 1.989e30  # Solar mass in kg
kpc_to_m = 3.086e19  # Kiloparsec to meters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
a0 = 1.2e-10  # MOND acceleration threshold (m/s^2)

# Functions to convert R.A. and Decl. to degrees
def hms_to_deg(h, m, s):
    return h * 15 + m / 4 + s / 240

def dms_to_deg(d, m, s):
    sign = 1 if d >= 0 else -1
    return sign * (abs(d) + m / 60 + s / 3600)


"""
#Loading the galaxy cluster
"""

galaxy_data = [
    {"mass": 8.29e7 * M_sun, "position": [-77, 0], "velocity": 2618, "ra_hms": (12, 43, 59), "dec_dms": (62, 20, 0)},  # D1
    {"mass": 2.75e8 * M_sun, "position": [0, 0], "velocity": 2624, "ra_hms": (12, 44, 12), "dec_dms": (62, 14, 51)},  # D2
    {"mass": 7.65e7 * M_sun, "position": [10.3, -10.3], "velocity": 2567, "ra_hms": (12, 44, 12), "dec_dms": (62, 10, 19)},  # D3
    {"mass": 1.47e7 * M_sun, "position": [20.6, -20.6], "velocity": 2550, "ra_hms": (12, 44, 20), "dec_dms": (62, 9, 58)},  # D4
    {"mass": 7.87e7 * M_sun, "position": [99, -77], "velocity": 2610, "ra_hms": (12, 44, 23), "dec_dms": (62, 3, 6)},  # D5
]

# Convert R.A. and Decl. to degrees and then to radians
for galaxy in galaxy_data:
    galaxy["ra_deg"] = hms_to_deg(*galaxy["ra_hms"])
    galaxy["dec_deg"] = dms_to_deg(*galaxy["dec_dms"])

# Reference velocity (D2's velocity) for relative calculations
reference_velocity = galaxy_data[1]["velocity"]

# Compute relative velocities and their x- and y-components
for galaxy in galaxy_data:
    # Convert line-of-sight velocity to relative velocity
    galaxy["velocity_relative"] = galaxy["velocity"] - reference_velocity
    
    # Convert R.A. and Decl. to radians
    ra_rad = np.radians(galaxy["ra_deg"])
    dec_rad = np.radians(galaxy["dec_deg"])
    
    # Compute x- and y-components of velocity
    galaxy["velocity_x"] = galaxy["velocity_relative"] * np.cos(dec_rad) * np.cos(ra_rad)
    galaxy["velocity_y"] = galaxy["velocity_relative"] * np.cos(dec_rad) * np.sin(ra_rad)

# Convert positions from kpc to meters, velocities from km/s to m/s
for galaxy in galaxy_data:
    galaxy["position"] = np.array(galaxy["position"]) * kpc_to_m
    galaxy["velocity"] = np.array([galaxy["velocity_x"], galaxy["velocity_y"]]) * 1e3  # Convert to m/s

# Create Galaxy objects
galaxies = [GD.Galaxy(g["mass"], g["position"], g["velocity"]) for g in galaxy_data]

# Display results
print("Results: x- and y-components of velocity (in m/s):")
for galaxy in galaxy_data:
    print(f"Galaxy: RA={galaxy['ra_deg']:.3f}°, Dec={galaxy['dec_deg']:.3f}°")
    print(f"  Relative Velocity: {galaxy['velocity_relative']} km/s")
    print(f"  Velocity x: {galaxy['velocity_x']:.2f} km/s")
    print(f"  Velocity y: {galaxy['velocity_y']:.2f} km/s")
"""

galaxy_data = [ 

    {"mass": 7.87e7 * M_sun, "position": [-77, 0], "velocity": [0, -33]},  # D1 

    {"mass": 2.75e8 * M_sun, "position": [0, 0], "velocity": [0, 0]},      # D2 

    {"mass": 7.65e7 * M_sun, "position": [10.3, -10.3], "velocity": [17, -30]},  # D3 

    {"mass": 1.47e7 * M_sun, "position": [20.6, -20.6], "velocity": [15, -27]},  # D4 

    {"mass": 7.87e7 * M_sun, "position": [99, -77], "velocity": [20, -10]},  # D5 

] 

  

# Convert positions from kpc to meters, velocities from km/s to m/s 

for galaxy in galaxy_data: 

    galaxy["position"] = np.array(galaxy["position"]) * kpc_to_m 

    galaxy["velocity"] = np.array(galaxy["velocity"]) * 1e3 

  

# Create Galaxy objects 

galaxies = [GD.Galaxy(g["mass"], g["position"], g["velocity"]) for g in galaxy_data] 
"""