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
        self.acceleration = np.zeros(3, dtype=np.float64)  # in m/s^2
    
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

# Interpolation function
def nu(y):
    return np.sqrt(0.5+0.5*np.sqrt(1+4/y**2))

def total_acceleration(galaxies, mode):
    num_galaxies = len(galaxies)
    positions = np.array([g.position for g in galaxies])
    masses = np.array([g.mass for g in galaxies])
    accelerations = np.zeros((num_galaxies, 3), dtype=np.float64)
    
    a0_inv = 1 / a0 if mode == "mond" else None
    
    for i, galaxy in enumerate(galaxies):
        r_vecs = positions - positions[i]  # Vector from galaxy `i` to others
        r_mags = np.linalg.norm(r_vecs, axis=1)  # Magnitude of r_vecs
        valid = r_mags > 0  # Avoid division by zero for self-interaction
        r_vecs = r_vecs[valid]
        r_mags = r_mags[valid]
        masses_valid = masses[valid]
        
        # Newtonian acceleration
        g = np.sum(G * masses_valid[:, None] * (r_vecs / r_mags[:, None]**3), axis=0)
        
        # Apply MOND modification if needed
        if mode == "mond":
            g_magnitude = np.linalg.norm(g)
            g *= nu(g_magnitude * a0_inv) if g_magnitude > 0 else 1
        
        accelerations[i] = g
    
    """
    Corrections for (angular) momentum preservation
    """
    total_mass = np.sum(masses)
    center_of_mass = np.sum(masses[:, None] * positions, axis=0) / total_mass
    total_force = np.sum(masses[:, None] * accelerations, axis=0)
    
    # T_1 and moment of inertia tensor
    T_1 = np.sum(
        masses[:, None] * np.cross(positions, accelerations), axis=0
    )
    
    inertia_tensor = np.zeros((3, 3), dtype=np.float64)

    for i in range(len(masses)):
        outer_product = np.outer(positions[i], positions[i])
        inertia_tensor += masses[i] * outer_product    
    
    T = T_1 - np.cross(center_of_mass, total_force)
    B = np.matmul(np.linalg.inv(np.trace(inertia_tensor) * np.eye(3) - inertia_tensor), T)
    
    A = total_force / total_mass
    
    # Adjust accelerations for MOND
    if mode == "mond":
        for i in range(num_galaxies):
            accelerations[i] -= (A + np.cross(B, (positions[i] - center_of_mass)))
    
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
    center_of_mass=np.array([0,0,0],dtype='float64')
    total_mass=0
    for galaxy in galaxies:
        kinetic_energy+=galaxy.kinetic_energy()
        center_of_mass+=galaxy.mass*galaxy.position
        total_mass+=galaxy.mass
    
    return center_of_mass/total_mass


