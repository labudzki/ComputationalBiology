#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 10:57:23 2025

@author: andrealabudzki
"""

import numpy as np
import matplotlib.pyplot as plt
import os


# Function defining line v6 = v1 + D
def func_v1_d(x, d):
    return x+d

# Commands to save plots to desktop
save_plots=True
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")

# Parameters
D = 0.7
vmax = 0.3

# Defining x values
x_lim = 1.1
x_vals = np.linspace(0, x_lim, 50)

#getting y = x+D values
y_vals = func_v1_d(x_vals, D)

    
#Plotting figure 
plt.figure()
plt.plot(x_vals, y_vals, label=f'v6 = v1 + D')


# vmax line
plt.axhline(y=vmax, color='black', linestyle='--', label=f'v6 = vmax')

# Filling in feasible region
plt.fill_between(x_vals[y_vals < vmax], y_vals[y_vals < vmax], vmax, color='orange', alpha=0.5)

plt.xlim(0, x_lim)
plt.ylim(0, x_lim)
plt.xlabel('v1')
plt.ylabel('v6')
plt.legend()
plt.grid(True)
plt.title(f'Accessible region, D = {D}, vmax = {vmax}')
if save_plots:
    plt.savefig(os.path.join(desktop_path, f"accessible_regions_D_{D}_vmax_{vmax}_plot.png"), dpi=300)
    plt.close()
else:
    plt.show()

