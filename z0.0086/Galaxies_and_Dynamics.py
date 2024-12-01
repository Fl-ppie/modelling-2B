# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:19:54 2024

@author: Fardeen
To do:
    
"""

import numpy as np

# Constants
M_sun = 1.989e30  # Solar mass in kg
kpc_to_m = 3.086e19  # Kiloparsec to meters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
a0 = 1.2e-10  # MOND acceleration threshold (m/s^2)


"""
Creating the Galaxy-class and dynamics on galaxies
"""
class Galaxy:
    def __init__(self, mass, position, velocity):
        self.mass = mass  # in kg
        self.position = np.array(position, dtype=np.float64)  # in meters
        self.velocity = np.array(velocity, dtype=np.float64)  # in meters per second
        self.acceleration = np.zeros(2, dtype=np.float64)  # in m/s^2
    
    def __str__(self):
        return f"mass: {self.mass}\nposition: {self.position}\nvelocity: {self.velocity}\nacceleration: {self.acceleration}\n"

# Dark matter halo parameters
M_dyn = 6.02e10 * M_sun  # Total dynamical mass
M_baryonic = 2.6e9 * M_sun  # Total baryonic mass
M_dark = M_dyn - M_baryonic  # Dark matter mass
r_virial = 154 * kpc_to_m / 2  # Half of the spatial span (virial radius)

"""
# NFW profile for dark matter halo
def nfw_enclosed_mass(r, M_vir=M_dark, r_vir=r_virial):
    c = 10  # Concentration parameter, typical for halos
    r_s = r_vir / c  # Scale radius
    rho_s = M_vir / (4 * np.pi * r_s**3 * (np.log(1 + c) - c / (1 + c)))
    enclosed_mass = 4 * np.pi * rho_s * r_s**3 * (np.log(1 + r / r_s) - r / (r + r_s))
    return enclosed_mass


# Acceleration due to dark matter
def dark_matter_acceleration(position):
    r = np.linalg.norm(position)
    if r == 0:
        return np.zeros(2)
    enclosed_mass = nfw_enclosed_mass(r)
    return -G * enclosed_mass / r**2 * position / r


# Newtonian acceleration between two masses
def newtonian_acceleration(m1, m2, r_vec):
    r_mag = np.linalg.norm(r_vec)
    if r_mag == 0:
        return np.zeros(2)
    return -G * m1 * m2 / r_mag**2 * r_vec / r_mag

# MOND interpolation function
def mond_interpolation(a_newton, a0=a0):
    x = a_newton / a0
    return x / np.sqrt(1 + x**2)

# Total acceleration calculation
def total_acceleration(galaxies, mode="newtonian"):
    accelerations = []
    for galaxy in galaxies:
        total_force = np.zeros(2)
        for other_galaxy in galaxies:
            if galaxy != other_galaxy:
                r_vec = other_galaxy.position - galaxy.position
                F_vec = newtonian_acceleration(galaxy.mass, other_galaxy.mass, r_vec)
                
                if mode == "mond":
                    a_newton = np.linalg.norm(F_vec / galaxy.mass)
                    mu = mond_interpolation(a_newton)
                    F_vec *= mu

                total_force += F_vec

        # Add dark matter contribution if enabled
        if mode == "dark_matter":
            total_force += dark_matter_acceleration(galaxy.position) * galaxy.mass

        accelerations.append(total_force / galaxy.mass)
    
    return accelerations
"""
# Dark matter acceleration
def dark_matter_acceleration(position):
    r = np.linalg.norm(position)
    if r == 0:
        return np.zeros(2)
    enclosed_mass = M_dark * (r / r_virial)**2  # Simplified density profile
    return -G * enclosed_mass / r**2 * position / r

# Total acceleration combining Newtonian, MOND, and dark matter
def total_acceleration(galaxies, mode="simple"):
    accelerations = []
    for galaxy in galaxies:
        total_force = np.zeros(2)
        for other_galaxy in galaxies:
            if galaxy != other_galaxy:
                r_vec = other_galaxy.position - galaxy.position
                r_mag = np.linalg.norm(r_vec)
                if r_mag == 0:
                    continue
                
                # Newtonian gravitational force
                F_newton = G * galaxy.mass * other_galaxy.mass / r_mag**2
                F_vec = F_newton * r_vec / r_mag
                
                # Apply MOND modification if enabled
                if mode != "newtonian":
                    a_newton = F_newton / galaxy.mass
                    mu = a_newton / (a_newton + a0)  # Simple interpolation
                    F_vec *= mu
                
                total_force += F_vec
        
        # Add dark matter contribution
        total_force += dark_matter_acceleration(galaxy.position) * galaxy.mass
        accelerations.append(total_force / galaxy.mass)
    return accelerations

"""
Estimation via the Runge-Kutta 4th Order Integration (RK4)
"""
def runge_kutta4(galaxies, dt, mode="simple"):
    accelerations = total_acceleration(galaxies, mode=mode)

    k1_v = [g.velocity for g in galaxies]
    k1_p = [g.position for g in galaxies]

    for i, galaxy in enumerate(galaxies):
        galaxy.velocity += accelerations[i] * dt / 2
        galaxy.position += galaxy.velocity * dt / 2

    k2_v = [g.velocity for g in galaxies]
    k2_p = [g.position for g in galaxies]

    for i, galaxy in enumerate(galaxies):
        galaxy.velocity = k1_v[i] + accelerations[i] * dt / 2
        galaxy.position = k1_p[i] + galaxy.velocity * dt / 2

    k3_v = [g.velocity for g in galaxies]
    k3_p = [g.position for g in galaxies]

    for i, galaxy in enumerate(galaxies):
        galaxy.velocity = k2_v[i] + accelerations[i] * dt
        galaxy.position = k2_p[i] + galaxy.velocity * dt

    k4_v = [g.velocity for g in galaxies]
    k4_p = [g.position for g in galaxies]

    for i, galaxy in enumerate(galaxies):
        galaxy.velocity = (k1_v[i] + 2 * k2_v[i] + 2 * k3_v[i] + k4_v[i]) / 6
        galaxy.position = (k1_p[i] + 2 * k2_p[i] + 2 * k3_p[i] + k4_p[i]) / 6

