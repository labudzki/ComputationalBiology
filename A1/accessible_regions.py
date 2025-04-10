#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 10:57:23 2025

@author: andrealabudzki
"""

import numpy as np
import matplotlib.pyplot as plt
import os

save_plots=True
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")

# Parameters
D_vals = [0.2, 0.5]
vmax = 0.8
x_vals = np.linspace(0, 1, 15)

for D in D_vals:
    y_vals = x_vals + D

    # Keep only points where both x and y are in [0,1]
    mask = (y_vals >= 0) & (y_vals <= 1)
    x = x_vals[mask]
    y = y_vals[mask]

    plt.figure()

    # Plot the line
    plt.plot(x, y, label=f'y = x + {D}')
    
    # Horizontal vmax line
    plt.axhline(y=vmax, color='black', linestyle='--', label=f'y = {vmax}')

    # Shade area between y and vmax where y < vmax
    mask_below = y < vmax
    plt.fill_between(x[mask_below], y[mask_below], vmax, color='orange', alpha=0.3)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.title(f'Shaded Region: y = x + {D} < {vmax}')
    if save_plots:
        plt.savefig(os.path.join(desktop_path, f"accessible_regions_D_{D}_plot.png"), dpi=300)
        plt.close()
    else:
        plt.show()

