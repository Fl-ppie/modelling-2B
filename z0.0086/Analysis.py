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
from copy import deepcopy

def plot_formatting(lims=[[None,None],[None,None]]):    
    plt.xlabel('X Position (m)')
    plt.ylabel('Y Position (m)')
    plt.legend()
    plt.title('Trajectories of Dwarf Galaxies')
    plt.grid(True)
    
    plt.xlim(lims[0])
    plt.ylim(lims[1])

if __name__=='__main__':
    filepath = 'position_data/'

    dt=1e13
    total_time=1e18

    # Run both models
    positions_newtonian, CoM_newtonian = Sim.run_simulation(deepcopy(Data.galaxies), mode="newtonian", dt=dt, total_time=total_time)
    
    # Save positions for later use in animation
    #np.save(filepath+'positions_newtonian', positions_newtonian)

    # Plot Newtonian positions
    plt.figure(figsize=(12, 6))
    for i in range(len(Data.galaxies)):
        plt.plot(positions_newtonian[i][:, 0], positions_newtonian[i][:, 1], label=f"Galaxy {i+1} Newtonian")

    plot_formatting()
    plt.show()
    
    plt.scatter(CoM_newtonian[:,0],CoM_newtonian[:,1])
    plt.title("Center of mass over time Newtonian")
    plt.show()
    
    #positions_dark_matter = run_simulation(galaxies, mode="dark matter")
    
    positions_mond, CoM_mond = Sim.run_simulation(deepcopy(Data.galaxies), mode="mond", dt=dt, total_time=total_time)
    
    # Save positions for later use in animation
    #np.save(filepath+'positions_mond', positions_mond)

    # Plot mond positions
    plt.figure(figsize=(12, 6))
    for i in range(len(Data.galaxies)):
        plt.plot(positions_mond[i][:, 0], positions_mond[i][:, 1], label=f"Galaxy {i+1} mond")

    plot_formatting()
    plt.show()
    
    plt.scatter(CoM_mond[:,0],CoM_mond[:,1])
    plt.title("Center of mass over time mond")
    plt.show()
    
    """
    Plotting the results
    """
    # Plot both force types in one
    plt.figure(figsize=(12, 6))
    for i in range(len(Data.galaxies)):
        plt.plot(positions_newtonian[i][:, 0], positions_newtonian[i][:, 1], label=f"Galaxy {i+1} Newtonian")
        plt.plot(positions_mond[i][:, 0], positions_mond[i][:, 1], label=f"Galaxy {i+1} MOND", linestyle="--")
        #plt.plot(positions_dark_matter[i][:, 0], positions_dark_matter[i][:, 1], label=f"Galaxy {i+1} Dark Matter", linestyle=":")
        
    plot_formatting()
    plt.show()