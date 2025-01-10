# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:19:54 2024

@author: Fardeen
To do:
    -correctie voor verschuiving massamiddelpunt
    
    
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
    
    def r(self):
        return np.linalg.norm(self.position)
    
    def theta(self):
        return np.atan2(self.position[1],self.position[0])
    
    def kinetic_energy(self):
        return 1/2*self.mass*np.linalg.norm(self.velocity)**2
    
    def angular_velocity(self):
        #(-v_x*y+v_y*x)/r^2
        return (-self.velocity[0]*self.position[1]+self.velocity[1]*self.position[0])/self.r()**2

# Dark matter halo parameters
M_dyn = 6.02e10 * M_sun  # Total dynamical mass
M_baryonic = 2.6e9 * M_sun  # Total baryonic mass
M_dark = M_dyn - M_baryonic  # Dark matter mass
r_virial = 154 * kpc_to_m / 2  # Half of the spatial span (virial radius)
# Dark matter acceleration
def dark_matter_acceleration(position):
    r = np.linalg.norm(position)
    if r == 0:
        return np.zeros(2)
    enclosed_mass = M_dark * (r / r_virial)**2  # Simplified density profile
    return -G * enclosed_mass / r**2 * position / r

# Interpolation function
def nu(x):
    return x/(1+x)

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
                if mode == "mond":
                    a_newton = F_newton / galaxy.mass
                    F_vec *= nu(a_newton/a0)
                
                total_force += F_vec
        # Add dark matter contribution
        #total_force += dark_matter_acceleration(galaxy.position) * galaxy.mass
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

    kinetic_energy=0
    center_of_mass=np.array([0,0],dtype='float64')
    total_mass=0
    for galaxy in galaxies:
        kinetic_energy+=galaxy.kinetic_energy()
        center_of_mass+=galaxy.mass*galaxy.position
        total_mass+=galaxy.mass
    
    return center_of_mass/total_mass
