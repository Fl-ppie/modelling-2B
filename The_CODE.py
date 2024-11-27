# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 20:20:00 2024

@author: Fardeen

TO DO:
    - Debuggen
    - Deep MOND module implementeren
    - Modulaire code
    - Data import van galaxies beter werkbaar
    - NFW profile for dark matter halo beter werkbaar
"""

"""
Preparations
"""
import numpy as np
import matplotlib.pyplot as plt

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
Creating the Galaxy-class and dynamics on galaxies
"""
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


"""
Loading the galaxy cluster
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
galaxies = [Galaxy(g["mass"], g["position"], g["velocity"]) for g in galaxy_data]

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
galaxies = [Galaxy(g["mass"], g["position"], g["velocity"]) for g in galaxy_data]
"""

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
        runge_kutta4(galaxies, dt, mode=mode)
        for i, galaxy in enumerate(galaxies):
            positions[i].append(galaxy.position.copy())
    for i in range(len(galaxies)):
        positions[i] = np.array(positions[i])
    return positions

# Run both models
positions_newtonian = run_simulation(galaxies, mode="newtonian")

plt.figure(figsize=(12, 6))
for i in range(len(galaxies)):
    plt.plot(positions_newtonian[i][:, 0], positions_newtonian[i][:, 1], label=f"Galaxy {i+1} Newtonian")

plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.title('Trajectories of Dwarf Galaxies')
plt.grid(True)
plt.axis('equal')
plt.show()

#positions_dark_matter = run_simulation(galaxies, mode="dark matter")
positions_mond = run_simulation(galaxies, mode="mond")

plt.figure(figsize=(12, 6))
for i in range(len(galaxies)):
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
for i in range(len(galaxies)):
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


plt.figure(figsize=(12, 6))
for i in range(len(galaxies)):
    plt.plot(positions_mond[i][:, 0], positions_mond[i][:, 1], label=f"Galaxy {i+1} MOND", linestyle="--")

plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.title('Trajectories of Dwarf Galaxies')
plt.grid(True)
plt.axis('equal')
plt.show()

"""
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