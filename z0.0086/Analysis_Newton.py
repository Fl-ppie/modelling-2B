# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 19:46:49 2025

@author: Fardeen
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

    dt=1e16
    total_time=1e21

    positions_newtonian, CoM_newtonian, properties_Newt = Sim.run_simulation(deepcopy(Data.galaxies), mode="newtonian", dt=dt, total_time=total_time)
    
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
    
    time_axis = np.arange(0,int(total_time), int(dt))/(364.25*24*3600*1e9)
    
    fig, axs = plt.subplots(3, 1, figsize=(8, 12))

    # Subplot 1: Total Kinetic Energy over Time
    axs[0].plot(time_axis, properties_Newt[0], '.')
    axs[0].set_xlabel("Time $t$ (billion year)")
    axs[0].set_ylabel("Total Kinetic Energy (J)")
    axs[0].set_title("Total Kinetic Energy over Time")

    # Subplot 2: Total Linear Momentum over Time
    axs[1].plot(time_axis, properties_Newt[1], '.')
    axs[1].set_xlabel("Time $t$ (billion year)")
    axs[1].set_ylabel("Total Linear Momentum (J)")
    axs[1].set_title("Total Linear Momentum over Time")

    # Subplot 3: Total Angular Momentum over Time
    axs[2].plot(time_axis, properties_Newt[2], '.')
    axs[2].set_xlabel("Time $t$ (billion year)")
    axs[2].set_ylabel("Total Angular Momentum (J)")
    axs[2].set_title("Total Angular Momentum over Time")

    # Adjust layout
    plt.tight_layout()
    plt.show()
    
    