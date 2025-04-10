#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 10:57:23 2025

@author: andrealabudzki
"""

import numpy as np
import matplotlib.pyplot as plt
import os


def func_v1_d(x, d):
    return x+d

save_plots=False
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")

# Parameters
# D_vals = [0.2, 0.5]
D = 1.2
vmax = 0.8
x_vals = np.linspace(0, 1.5, 20)

#y = x+D
y_vals = func_v1_d(x_vals, D)

    

plt.figure()
plt.plot(x_vals, y_vals, label=f'v6 = v1 + D')
plt.plot(x_vals, x_vals, label=f'v6 = v1')

# vmax line
plt.axhline(y=vmax, color='black', linestyle='--', label=f'v6 = vmax')

plt.axvline(x=D, color='black', label=f'x=D')

plt.fill_between(x_vals[y_vals < vmax], y_vals[y_vals < vmax], vmax, color='orange', alpha=0.3)

plt.xlim(0, 1.5)
plt.ylim(0, 1.5)
plt.xlabel('v1')
plt.ylabel('v6')
plt.legend()
plt.grid(True)
plt.title(f'Accessible region, D = {D}')
if save_plots:
    plt.savefig(os.path.join(desktop_path, f"accessible_regions_D_{D}_plot.png"), dpi=300)
    plt.close()
else:
    plt.show()



# for D in D_vals:
#     y_vals1 = x_vals + D
#     y_vals2 = x_vals 

#     # Keep only points where both x and y are in [0,1]
#     mask = (y_vals1 >= 0) & (y_vals1 <= 1)
#     x = x_vals[mask]
#     y1 = y_vals1[mask]
#     y2 = y_vals2[mask]

#     plt.figure()
#     plt.plot(x, y1, label=f'v6 = v1 + D')
#     plt.plot(x, y2, label=f'v6 = v1')
    
#     # vmax line
#     plt.axhline(y=vmax, color='black', linestyle='--', label=f'v6 = vmax')

#     # Shade area between y and vmax where y < vmax
#     mask_below = y1 < vmax
#     plt.fill_between(x[mask_below], y1[mask_below], vmax, color='orange', alpha=0.3)

#     plt.xlim(0, 1)
#     plt.ylim(0, 1)
#     plt.xlabel('v1')
#     plt.ylabel('v6')
#     plt.legend()
#     plt.grid(True)
#     plt.title(f'Accessible region')
#     if save_plots:
#         plt.savefig(os.path.join(desktop_path, f"accessible_regions_D_{D}_plot.png"), dpi=300)
#         plt.close()
#     else:
#         plt.show()

