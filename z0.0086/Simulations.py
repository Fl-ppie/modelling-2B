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
def run_simulation(galaxies, mode, dt, total_time):
    steps = int(total_time / dt)
    positions = [[] for i in range(len(galaxies))]
    CoM = []
    kin_energy = []
    lin_momentum = []
    angular_momentum = []
    
    
    for step in range(steps):
        if step % 100 == 0:
            print(f"Step {step}/{steps}")
            
        CoM.append(GD.runge_kutta4(galaxies, dt, mode=mode))
        total_kin = 0.0
        total_lin_moment = 0.0
        total_ang_moment = 0.0
        
        for i, galaxy in enumerate(galaxies):
            positions[i].append(np.copy(galaxy.position))
            total_kin += galaxy.kinetic_energy()
            total_lin_moment += galaxy.mass*np.copy(galaxy.velocity)
            total_ang_moment += galaxy.mass*np.cross(np.copy(galaxy.position), np.copy(galaxy.velocity))
        
        kin_energy.append(total_kin)
        lin_momentum.append(np.linalg.norm(total_lin_moment))
        angular_momentum.append(np.linalg.norm(total_ang_moment))
        
    properties = [kin_energy, lin_momentum, angular_momentum]
    return np.array(positions), np.array(CoM), properties
