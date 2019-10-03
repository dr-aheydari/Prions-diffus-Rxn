#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Created on Fri Sep 20 09:31:58 2019
    
    @author: aliheydari
    """

## PLOT THE Mass of A against B

import numpy as np
# importing the required module
import matplotlib.pyplot as plt



path_psi = "/Users/aliheydari/Documents/Prions Project Data/total_psi_data_diffOnly/"
path_zeta = "/Users/aliheydari/Documents/Prions Project Data/total_zeta_data_diffOnly/"

mass_A = np.loadtxt(path_psi + "A_SoftDist_DStays2.txt");
mass_B = np.loadtxt(path_zeta + "B_SoftDist_DStays2.txt");

# only get the values not the time-steps

mass_A = mass_A[:,1];
mass_B = mass_B[:,1];

if len(mass_A) == len(mass_B):
    
    x = np.arange(len(mass_A));
else :
    print("I fucked up!");


# plotting the points
SpeciesA = plt.plot(x, mass_A, color='green',label='A', linestyle='dashed', linewidth = 3
                    ,marker='o', markerfacecolor='blue', markersize=5)

SpeciesB = plt.plot(x,mass_B, color='red',label='B', linewidth = 3
                    ,marker='o', markerfacecolor='yellow', markersize=5)

#
# naming the x axis
plt.xlabel('Time',fontsize=20)
# naming the y axis
plt.ylabel('Mass',fontsize=20)

# giving a title to my graph
plt.title('Mass of Species A and B vs. Time',fontsize=14)
# plotting the points
plt.legend(('Species A', 'Species B'))


plt.savefig('masses_DiffReaction.png',dpi = 1080)


plt.show()


########### TOTAL MASS ############

total_mass = mass_A + mass_B;

# to set the precision for a vector

float_formatter = lambda x: "%.6f" % x

# change the precision to undermine numerical error, if you want
for i in range(len(total_mass)):
    total_mass[i] = float_formatter(total_mass[i])


# plot
totMass = plt.plot(x, total_mass, color='purple',label='A', linestyle='dashed', linewidth = 5
                   ,marker='o', markerfacecolor='blue', markersize=5)

plt.xlabel('Time',fontsize=16)
# naming the y axis
plt.ylabel('Total Mass ',fontsize=16)

# giving a title to my graph
plt.title('Total Mass (Order of 1e-6) ',fontsize=18)

plt.savefig('TotMass_DiffReaction.png',dpi = 1080)


plt.show() 
