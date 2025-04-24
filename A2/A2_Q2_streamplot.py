import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy.optimize as sp

"""
Forked and modified (butchered) code from viterbi.py
New functions:
- gene_regulation_model_3
Modified functions:
- 'vector_field'
- 'plot_phase_plane'
"""

def gene_regulation_model_3(y0, t, alpha, beta, gamma, delta): 
    """
    Gene regulation model for ...

    Parameters: description listed in the main() function. 

    Returns: the derivatives of the concentrations of metabolite and protein.
    """
    x, y = y0

    dx_dt = alpha * x - beta*x*y
    dy_dt = -gamma * y + delta * x * y

    return dx_dt, dy_dt


def solve_gene_regulation_model(gene_regulation_model, initial_conditions, t, params):
    """
    Solves the specified gene regulation model using the odeint solver.
    """
    solution = odeint(gene_regulation_model, initial_conditions, t, args=(
        params['alpha'], params['beta'], params['gamma'], params['delta']
    ))
    return solution


def vector_field(model_type, grid, initial_conditions, params):
    """
    Computes the vector field for the phase plane plot.

    Parameters: 
    model_type: function
        The gene regulation model function to be used.
    grid: np.ndarray
        The grid of values for the phase plane.
    initial_conditions: list
        Initial conditions for the model.
    params: dict
        Parameters for the gene regulation model.

    Returns:
    U: np.ndarray
        The x values of the vector field for the phase plot.
    V: np.ndarray
        The y values of the vector field for the phase plot.
    """
    U = np.zeros((len(grid), len(grid)))
    V = np.zeros((len(grid), len(grid)))
    
    for i, x in enumerate(grid):
        for j, y in enumerate(grid):
            y0 = [x, y]
            dydt = model_type(y0, 0, **params)
            U[j, i] = dydt[0]  # dx/dt
            V[j, i] = dydt[1]  # dy/dt
    
    print("U max:", np.max(np.abs(U)))
    print("V max:", np.max(np.abs(V)))
    return U, V


def plot_phase_plane(model_type, initial_conditions, var_indices, grid, params, solution=None, nullclines=False, equilibria=False, save_path=None):
    """
    Plots the phase plane for the gene regulation model.

    Parameters:
    model_type: function
        The gene regulation model function to be used.
    initial_conditions: list
        Initial conditions for the model.
    var_indices: list
        Indices of the variables to be plotted on the x and y axes.
    grid: np.ndarray
        The grid of values for the phase plane.
    params: dict
        Parameters for the gene regulation model.
    solution: np.ndarray
        Solution array for plotting the trajectory, solved by odeint solver.
    save_path: str
        Path to save the plot. If None, the plot will be shown.
    """
    X, Y = np.meshgrid(grid, grid)
    U, V = vector_field(model_type, grid, initial_conditions, params)
    speed = np.sqrt(U**2 + V**2) # Magnitude of the vector field

    plt.figure(figsize=(8, 6))
    strm = plt.streamplot(X, Y, U, V, color=speed, cmap='magma', density=1.0, linewidth=1)
    plt.colorbar(strm.lines, label=r'Speed (|$dc/dt$|)')

    plt.xlabel(r'Metabolite X ($x$)')
    plt.ylabel(r'Enzyme Y ($y$)')
    plt.title(r'Phase Plane: $x$ vs $y$ (colored by speed)')
    plt.grid(True)

    if solution is not None:
        plt.plot(solution[:, var_indices[0]], solution[:, var_indices[1]], color='black', lw=2, label='Trajectory')
        plt.legend()

    if nullclines:
        # x-nullcline
        plt.axvline(x=0, color='red', lw=2, linestyle=":", label='x-nullcline(s)')
        plt.axhline(y=params["alpha"]/params["beta"], color='red', lw=2, linestyle=":")

        # y-nullcline
        plt.axvline(x=params["gamma"]/params["delta"], color='blue', lw=2, linestyle=":", label='y-nullcline')
        plt.axhline(y=0, color='blue', lw=2, linestyle=":")        
        plt.legend()

    if equilibria:
        pnts = [(0,0), (params["gamma"]/params["delta"], params["alpha"]/params["beta"])]
        plt.scatter(*zip(*pnts), color='black', label="Equilibrium points")
        plt.legend(loc='upper right')

    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()

initial_conditions = [1, 0.5]  # Initial concentrations for x, y
t = np.linspace(0, 100, 1000)  # small steps for smooth limit-cycle

params = {
    'alpha': 2, # metabolite import rate
    'beta': 1.1, # complex formation rate
    'gamma': 1, # enzyme degradation rate
    'delta': 0.9 # catalysis rate
}

solution = solve_gene_regulation_model(
    gene_regulation_model_3,
    initial_conditions,
    t,
    params
)

# Plot phase plane for model 1
var_indices = [0,1]  # x vs y
grid_vals = np.linspace(-2.5, 10, 20) # show saddle-point at (0,0)

plot_phase_plane(
    gene_regulation_model_3,
    var_indices=var_indices,
    initial_conditions=initial_conditions,
    grid=grid_vals,
    params=params,
    nullclines=True,
    equilibria=True,
    solution=solution
)