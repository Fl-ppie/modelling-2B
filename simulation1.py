import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11 # Gravitational constant (m^3 kg^-1 s^-2)
a_0 = 1.2e-10   # Milgrom's constant (m/s^2)
M_0 = 1.988e30  # Solar mass (kg)
ly = 9.416e15   # Light year (m)

path = 'images/frames/'

# Define the Galaxy class
class Galaxy:
    def __init__(self, mass, position, velocity):
        self.mass = mass  # in kg
        self.position = np.array(position, dtype=np.float64)  # in meters
        self.velocity = np.array(velocity, dtype=np.float64)  # in meters per second
        self.acceleration = np.zeros(2, dtype=np.float64)  # in m/s^2

# Gravitational force between two galaxies
def gravitational_force(galaxy1, galaxy2):
    r_vec = galaxy2.position - galaxy1.position
    r_mag = np.linalg.norm(r_vec)
    if r_mag == 0:
        return np.zeros(2)
    
    #force_magnitude = G * galaxy1.mass * galaxy2.mass / r_mag**2 #NEWTON
    force_magnitude = 2 / 3 * (G*a_0)**(1/2) * ((galaxy1.mass + galaxy2.mass)**(3/2) - galaxy1.mass**(3/2) - galaxy2.mass**(3/2)) / r_mag #MOND
    
    force_vector = force_magnitude * r_vec / r_mag  # Force vector pointing towards the other galaxy
    return force_vector

# Calculate the acceleration of each galaxy based on gravitational forces
def calculate_accelerations(galaxies):
    accelerations = []
    for galaxy in galaxies:
        total_force = np.zeros(2)
        for other_galaxy in galaxies:
            if galaxy != other_galaxy:
                force = gravitational_force(galaxy, other_galaxy)
                total_force += force
        acceleration = total_force / galaxy.mass  # a = F/m
        accelerations.append(acceleration)
    return accelerations

# Runge-Kutta 4th Order Integration (RK4)
def runge_kutta4(galaxies, dt):
    # Calculate initial accelerations
    accelerations = calculate_accelerations(galaxies)
    
    # Calculate k1, k2, k3, k4 for positions and velocities
    k1_v = [g.velocity for g in galaxies]
    k1_p = [g.position for g in galaxies]
    
    # Calculate intermediate velocities and positions
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
    
    # Update positions and velocities using the Runge-Kutta method
    for i, galaxy in enumerate(galaxies):
        galaxy.velocity = (k1_v[i] + 2 * k2_v[i] + 2 * k3_v[i] + k4_v[i]) / 6
        galaxy.position = (k1_p[i] + 2 * k2_p[i] + 2 * k3_p[i] + k4_p[i]) / 6

def plotting():
    plt.xlabel('X position (meters)')
    plt.ylabel('Y position (meters)')
    plt.title('Galaxies in Local Cluster Simulation (MOND)')
    plt.legend()
    plt.grid(True)
    plt.xlim([-1e22,1e22])
    plt.ylim([-1e22,1e22])

# Initialize the galaxies (mass in kg, position in meters, velocity in m/s)
num_galaxies = 2  # Number of galaxies in the cluster
galaxies = []

# Example: Initialize galaxies in a more stable configuration (closer, slower)
np.random.seed(42)  # Set seed for reproducibility
for i in range(num_galaxies):
    mass = np.random.uniform(1e9*M_0, 1e9*M_0)
    position = np.random.uniform(-1e6*ly, 1e6*ly, size=2)
    velocity = np.random.uniform(-1e4, 1e4, size=2)
    galaxies.append(Galaxy(mass, position, velocity))

# Time parameters
dt = 1e14  # Time step
total_time = 1e19  # Total time
steps = int(total_time / dt)
printstep = 1e2

# List for storing positions for plotting later
positions = {i: [] for i in range(len(galaxies))}

# Run the simulation for the specified time
for step in range(steps):
    runge_kutta4(galaxies, dt)
    for i, galaxy in enumerate(galaxies):
        positions[i].append(galaxy.position.copy())
        
    # Save an image every so often for an animation
    if step%printstep==0 and step != 0:
        temp = {i: [] for i in range(len(galaxies))}
        for i in range(len(galaxies)):
            temp[i] = np.array(positions[i])
    
        plt.figure(figsize=(8, 8))
        for i in range(len(galaxies)):
            plt.plot(temp[i][-1000:, 0], temp[i][-1000:, 1], label=f"Galaxy {i+1}")
            plt.scatter(temp[i][-1,0], temp[i][-1,1])
        
        plotting()
        temp = str(int(step/printstep))
        name = (len(str(int(steps/printstep)))-len(temp))*'0'+temp
        plt.savefig(path+f'{name}.png')
        plt.close()
        print(step)

# Convert positions to numpy arrays for easy plotting
for i in range(len(galaxies)):
    positions[i] = np.array(positions[i])

# Plot the trajectory of the galaxies
plt.figure(figsize=(8, 8))
for i in range(len(galaxies)):
    plt.plot(positions[i][:, 0], positions[i][:, 1], label=f"Galaxy {i+1}")
    plt.scatter(positions[i][-1,0], positions[i][-1,1])

plotting()
plt.savefig(path+f'{int((step+1)/printstep)}.png')
plt.show()