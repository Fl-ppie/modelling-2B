# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:19:55 2024

@author: Fardeen

To Do:
    - 
"""
import Cluster_Data as Data
import Simulations as Sim
import numpy as np
import matplotlib.pyplot as plt


dt=1e13
total_time=1e18

# Run both models
positions_newtonian = Sim.run_simulation(Data.galaxies, mode="newtonian", dt=dt, total_time=total_time)

plt.figure(figsize=(12, 6))
for i in range(len(Data.galaxies)):
    plt.plot(positions_newtonian[i][:, 0], positions_newtonian[i][:, 1], label=f"Galaxy {i+1} Newtonian")

plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.title('Trajectories of Dwarf Galaxies')
plt.grid(True)
plt.axis('equal')
plt.show()

#positions_dark_matter = run_simulation(galaxies, mode="dark matter")
positions_mond = Sim.run_simulation(Data.galaxies, mode="mond", dt=dt, total_time=total_time)

plt.figure(figsize=(12, 6))
for i in range(len(Data.galaxies)):
    plt.plot(positions_mond[i][:, 0], positions_mond[i][:, 1], label=f"Galaxy {i+1} MOND", linestyle="--")

plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.title('Trajectories of Dwarf Galaxies')
plt.grid(True)
plt.axis('equal')
plt.show()



"""
Plotting the results
"""
plt.figure(figsize=(12, 6))
for i in range(len(Data.galaxies)):
    plt.plot(positions_newtonian[i][:, 0], positions_newtonian[i][:, 1], label=f"Galaxy {i+1} Newtonian")
    plt.plot(positions_mond[i][:, 0], positions_mond[i][:, 1], label=f"Galaxy {i+1} MOND", linestyle="--")
    #plt.plot(positions_dark_matter[i][:, 0], positions_dark_matter[i][:, 1], label=f"Galaxy {i+1} Dark Matter", linestyle=":")

plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.title('Trajectories of Dwarf Galaxies')
plt.grid(True)
plt.axis('equal')
plt.show()

"""
plt.figure(figsize=(12, 6))
for i in range(len(Data.galaxies)):
    plt.plot(positions_mond[i][:, 0], positions_mond[i][:, 1], label=f"Galaxy {i+1} MOND", linestyle="--")

plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.title('Trajectories of Dwarf Galaxies')
plt.grid(True)
plt.axis('equal')
plt.show()


plt.figure(figsize=(12, 6))
for i in range(len(galaxies)):
    plt.plot(positions_dark_matter[i][:, 0], positions_dark_matter[i][:, 1], label=f"Galaxy {i+1} Dark Matter", linestyle=":")

plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.title('Trajectories of Dwarf Galaxies')
plt.grid(True)
plt.axis('equal')
plt.show()
"""

"""
SOME OLD DATA - IGNORE

class Galaxy:
    def __init__(self, mass, position, velocity):
        self.mass = mass  # in kg
        self.position = np.array(position, dtype=np.float64)  # in meters
        self.velocity = np.array(velocity, dtype=np.float64)  # in meters per second
        self.acceleration = np.zeros(2, dtype=np.float64)  # in m/s^2

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
galaxies = [Galaxy(g["mass"], g["position"], g["velocity"]) for g in galaxy_data]
"""