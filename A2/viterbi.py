#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 15:31:20 2025

@author: andrealabudzki
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy.optimize as sp
import os


def viterbi_algorithm(obs, states, start_probs, trans_probs, emit_probs):
    """
    Viterbi algorithm for finding the most probable sequence of hidden states
    given a sequence of nucleotides.

    Parameters:
    obs: list
        Observations (consists of list of nucleotides in give string, ex. ACAGT).
    states: list
        Hidden states (ex. Exon, Intron).
    start_probs: dict
        Initial probabilities of each state (Pi).
    trans_probs: dict
        Transition probabilities between states.
    emit_probs: dict
        Emission probabilities of observations given states.
    
    Returns:
    final_prob: float
        Probability of the most probable path.
    best_path: list
        The most probable path of hidden states given the final state,
        based on the maximum final probability.
    """
    dp = [{}] # dp[t][state] = maximum probability of state at time t
    backtrack = {} # Path to the current state

    for state in states:
        dp[0][state] = start_probs[state] * emit_probs[state][obs[0]] # Calculates the initial probabilities
        backtrack[state] = [state]

    for t in range(1, len(obs)):
        dp.append({})
        new_backtrack = {}
        
        # Calculating the remaining probabilities for remaining observations
        for current_state in states:
            (max_prob, best_prev_state) = max(
                (
                    dp[t - 1][prev_state] * trans_probs[prev_state][current_state] * emit_probs[current_state][obs[t]],
                    prev_state 
                )
                for prev_state in states
            )

            dp[t][current_state] = max_prob
            new_backtrack[current_state] = backtrack[best_prev_state] + [current_state] # Keeps track of the best path to the current state (Psi)
        backtrack = new_backtrack

    final_state, final_prob = max(dp[-1].items(), key=lambda item: item[1])
    best_path = backtrack[final_state] # Best path to final state using backtracking 

    return final_prob, best_path

def gene_regulation_model_1(y, t, m_a, m_b, n_a, n_b, theta_a, theta_b, gamma_a, gamma_b, k_a, k_b, delta_a, delta_b): 
    """
    Gene regulation model for mRNA and protein concentrations for route I. 

    Parameters: description listed in the main() function. 

    Returns: the derivatives of the concentrations of mRNA and proteins.
    """
    r_a, r_b, p_a, p_b = y

    dr_a_dt = m_a * (p_b**n_b)/(p_b**n_b + theta_b**n_b) - gamma_a * r_a
    dr_b_dt = m_b * (theta_a**n_a)/(p_a**n_a + theta_a**n_a) - gamma_b * r_b

    dp_a_dt = k_a * r_a - delta_a * p_a
    dp_b_dt = k_b * r_b - delta_b * p_b
    return dr_a_dt, dr_b_dt, dp_a_dt, dp_b_dt

def gene_regulation_model_2(y, t, dt, m_a, m_b, n_a, n_b, theta_a, theta_b, gamma_a, gamma_b, k_a, k_b, delta_a, delta_b, sigma_1a, sigma_2a, sigma_1b, sigma_2b, noise=True, dt=1e-2):
    """
    Stochastic gene regulation model using SDEVelo formalism (Route II).
    y = [u_a, u_b, s_a, s_b, p_a, p_b] (unspliced, spliced mRNAs, proteins)
    """

    u_a, u_b, s_a, s_b, p_a, p_b = y
    du_a_dt = m_a * (p_b**n_b)/(p_b**n_b + theta_b**n_b) - gamma_a * u_a
    du_b_dt = m_b * (theta_a**n_a)/(p_a**n_a + theta_a**n_a) - gamma_b * u_b
    ds_a_dt = gamma_a * u_a - k_a * s_a
    ds_b_dt = gamma_b * u_b - k_b * s_b
    dp_a_dt = k_a * s_a - delta_a * p_a
    dp_b_dt = k_b * s_b - delta_b * p_b
    if noise:
        du_a_dt += sigma_1a * np.random.normal(0, np.sqrt(dt))
        du_b_dt += sigma_1b * np.random.normal(0, np.sqrt(dt))
        ds_a_dt += sigma_2a * np.random.normal(0, np.sqrt(dt))
        ds_b_dt += sigma_2b * np.random.normal(0, np.sqrt(dt))

    return du_a_dt, du_b_dt, ds_a_dt, ds_b_dt, dp_a_dt, dp_b_dt

def solve_gene_regulation_model(gene_regulation_model, initial_conditions, t, params):
    """
    Solves the specified gene regulation model using the odeint solver.
    """
    solution = odeint(gene_regulation_odel, initial_conditions, t, args=(
        params['m_a'], params['m_b'], params['n_a'], params['n_b'], 
        params['theta_a'], params['theta_b'], params['gamma_a'], params['gamma_b'], 
        params['k_a'], params['k_b'], params['delta_a'], params['delta_b']
    ))
    return solution

def plot_mRNA_time_evolution(t, solution, save_path=None):
    """
    Plots the results of the gene regulation model.
    """
    r_A = solution[:, 0]
    r_B = solution[:, 1]
    plt.figure(figsize=(10, 6))
    plt.plot(t, r_A, label='mRNA A', color='blue')
    plt.plot(t, r_B, label='mRNA B', color='red')
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (M)')
    plt.title('Gene Regulation Model')
    plt.legend()
    plt.grid()
    
    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()

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
    
    r_A, r_B = initial_conditions[0], initial_conditions[1]

    for i, p_A in enumerate(grid):
        for j, p_B in enumerate(grid):
            y0 = [r_A, r_B, p_A, p_B]
            dydt = model_type(y0, 0, **params)
            U[j, i] = dydt[2]  # dp_A/dt
            V[j, i] = dydt[3]  # dp_B/dt

    
    print("U max:", np.max(np.abs(U)))
    print("V max:", np.max(np.abs(V)))
    return U, V


def plot_phase_plane(model_type, initial_conditions, var_indices, grid, params, solution=None, save_path=None):
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
    plt.colorbar(strm.lines, label=r'Speed (|$dp/dt$|)')

    plt.xlabel(r'Protein A ($p_a$)')
    plt.ylabel(r'Protein B ($p_b$)')
    plt.title(r'Phase Plane: $p_a$ vs $p_b$ (colored by speed)')
    plt.grid(True)

    if solution is not None:
        plt.plot(solution[:, var_indices[0]], solution[:, var_indices[1]], color='black', lw=2, label='Trajectory')
        plt.legend()

    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()

def main():
    """
    Main function to run the program for the Viterbi algorithm and gene regulation models. 
    """
    # Setting up libraries and arrays to represent the values in the tables given in the task
    states = ['Exon', 'Intron']
    observations = ['A', 'G', 'C', 'G', 'C']

    start_probabilities = {
        'Exon': 0.5,
        'Intron': 0.5
    }

    transition_probabilities = {
        'Exon': {'Exon': 0.9, 'Intron': 0.1},
        'Intron': {'Exon': 0.2, 'Intron': 0.8}
    }

    emission_probabilities = {
        'Exon': {'A': 0.25, 'U': 0.25, 'G': 0.25, 'C': 0.25},
        'Intron': {'A': 0.4, 'U': 0.4, 'G': 0.05, 'C': 0.15}
    }

    # Run Veterbi algorithm
    probability, state_path = viterbi_algorithm(
        observations,
        states,
        start_probabilities,
        transition_probabilities,
        emission_probabilities
    )

    print("Most probable state path:", state_path)
    print("Probability of the path:", probability)

    initial_conditions = [0.8, 0.8, 0.8, 0.8]  # Initial concentrations for r_A, r_B, p_A, p_B
    t = np.linspace(0, 100, 100)  

    det_params = {
        'm_a': 2.35, 'm_b': 2.35,  # Max transcription rates
        'n_a': 3, 'n_b': 3,        # Hill coefficients
        'theta_a': 0.21, 'theta_b': 0.21,  # Binding thresholds
        'gamma_a': 1.0, 'gamma_b': 1.0,    # mRNA degradation rates
        'k_a': 1.0, 'k_b': 1.0,    # Translation rates
        'delta_a': 1.0, 'delta_b': 1.0  # Protein degradation rates
    }

    stoch_params = {
        'a_A': 1.0, 'a_B': 0.25,
        'b_A': 0.0005, 'b_B': 0.0005,
        'c_A': 2.0, 'c_B': 0.5,
        'beta_A': 2.35, 'beta_B': 2.35,
        'gamma_A': 1.0, 'gamma_B': 1.0,
        'n_A': 3, 'n_B': 3,
        'theta_A': 0.21, 'theta_B': 0.21,
        'k_PA': 1.0, 'k_PB': 1.0,
        'delta_PA': 1.0, 'delta_PB': 1.0,
        'sigma_1A': 0.05, 'sigma_2A': 0.05,
        'sigma_1B': 0.05, 'sigma_2B': 0.05,
    }

    GRN1_solution = solve_gene_regulation_model(
        gene_regulation_model_1,
        initial_conditions,
        t,
        det_params
    )

    plot_mRNA_time_evolution(t, GRN1_solution)
    
    # Plot phase plane for model 1
    var_indices = [2, 3]  # p_A vs p_B
    grid_vals = np.linspace(0, 3, 20)

    plot_phase_plane(
        gene_regulation_model_1,
        var_indices=var_indices,
        initial_conditions=initial_conditions,
        grid=grid_vals,
        params=det_params,
        solution=GRN1_solution
    )





if __name__ == "__main__":
    main()
