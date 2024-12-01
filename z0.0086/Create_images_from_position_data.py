# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 14:45:51 2024

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt

from Analysis import plot_formatting

def get_images(positions, savepath, mode, total_images=500):
    length = positions.shape[1]
    interval = int(length/total_images)
    
    tail_length = interval*10
    
    maximum = 0
    minimum = 0
    
    # Determine maximum and minimum values reached
    for i in range(len(positions)):
        for j in range(len(positions[i])):
            for k in range(len(positions[i,j])):
                if positions[i,j,k] > maximum:
                    maximum = positions[i,j,k]
                elif positions[i,j,k] < minimum:
                    minimum = positions[i,j,k]
    
    size = 0
    if maximum>-minimum:
        size = maximum
    if minimum<-maximum:
        size = -minimum
    
    # Add a margin
    size *= 1.05
    
    # One loop per frame
    for i in range(0, total_images):
        # Create 'tail'
        pos = positions[:,i*interval:i*interval+tail_length]
        number = (len(str(total_images))-len(str(i)))*'0'+str(i)
        
        # Plot 'tail' with dot at the head
        plt.figure(figsize=(8,8))
        for j in range(len(pos)):
            plt.plot(pos[j,:,0],pos[j,:,1], label=f"Galaxy {j+1} {mode}")
            plt.scatter(pos[j,-1,0],pos[j,-1,1], s=15)
        
        plot_formatting([[-size,size],[-size,size]])
        plt.savefig(savepath+number+'.png')
        plt.close()
        
        print(f"Created image {i+1}/{total_images}")
    
    plt.figure(figsize=(8,8))
    
    for i in range(len(positions)):
        plt.plot(positions[i,:,0],positions[i,:,1], label=f"Galaxy {i+1} {mode}")
        plt.scatter(positions[i,-1,0],positions[i,-1,1])
    
    plot_formatting([[-size,size],[-size,size]])
    plt.savefig(savepath+str(total_images)+'.png')
    plt.close()
        
if __name__=='__main__':
    
    # Location of files
    filepath = 'position_data/'
    path_newton = 'positions_newtonian.npy'
    path_mond = 'positions_mond.npy'
    
    # Load data
    pos_newton = np.load(filepath+path_newton)
    pos_mond = np.load(filepath+path_mond)
    
    # Create the images
    get_images(pos_newton, filepath+'images_newtonian/', 'Newtonian')
    get_images(pos_mond, filepath+'images_mond/', 'MOND')