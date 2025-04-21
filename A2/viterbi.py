#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stochastic & Deterministic Gene Regulation Model + Viterbi HMM

Created on Apr 17, 2025
@author: andrealabudzki
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def viterbi_algorithm(obs, states, start_probs, trans_probs, emit_probs):
    """
    Viterbi algorithm for finding the most probable sequence of hidden states 
    given a nucleotide sequence.
    """
    dp = [{}]
    backtrack = {}

    for state in states:
        dp[0][state] = start_probs[state] * emit_probs[state][obs[0]]
        backtrack[state] = [state]

    for t in range(1, len(obs)):
        dp.append({})
        new_backtrack = {}
        for curr_state in states:
            (max_prob, prev_best) = max(
                (
                    dp[t - 1][prev_state] * trans_probs[prev_state][curr_state] * emit_probs[curr_state][obs[t]],
                    prev_state
                ) for prev_state in states
            )
            dp[t][curr_state] = max_prob
            new_backtrack[curr_state] = backtrack[prev_best] + [curr_state]
        backtrack = new_backtrack

    final_state, final_prob = max(dp[-1].items(), key=lambda x: x[1])
    return final_prob, backtrack[final_state]


def gene_regulation_det(y, t, m_a, m_b, n_a, n_b, theta_a, theta_b, gamma_a, gamma_b, k_a, k_b, delta_a, delta_b):
    """
    Deterministic ODE model of gene regulation with Protein A inhibiting transcription of Gene B, 
    and Protein B promoting transcription of Gene A.
    """
    r_a, r_b, p_a, p_b = y

    # Transcription and splicing
    dr_a_dt = m_a * (p_b**n_b) / (p_b**n_b + theta_b**n_b) - gamma_a * r_a
    dr_b_dt = m_b * (theta_a**n_a) / (p_a**n_a + theta_a**n_a) - gamma_b * r_b

    # Translation
    dp_a_dt = k_a * r_a - delta_a * p_a
    dp_b_dt = k_b * r_b - delta_b * p_b

    return dr_a_dt, dr_b_dt, dp_a_dt, dp_b_dt


def solve_gene_regulation_det(model, initial_conditions, t, params):
    """
    Solve the deterministic ODE model using scipy's odeint.
    """
    return odeint(model, initial_conditions, t, args=tuple(params.values()))


def plot_mRNA_time_evolution_det(t, solution, save_path=None):
    """
    Plot the deterministic time evolution of mRNA concentrations for A and B (Route I). 
    """
    r_a, r_b = solution[:, 0], solution[:, 1]
    plt.figure(figsize=(10, 6))
    plt.plot(t, r_a, label='mRNA A', color='blue')
    plt.plot(t, r_b, label='mRNA B', color='red')
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration')
    plt.title('Deterministic Gene Regulation')
    plt.legend()
    plt.grid()
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()


def vector_field_det(model_type, grid, initial_conditions, params):
    """
    Calculate the vector field for the deterministic model (Route I). 
    This function computes the derivatives of the protein concentrations
    at each point in the grid to plot the phase plot. 
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


def plot_phase_plane_det(model_type, initial_conditions, var_indices, grid, params, solution=None, save_path=None):
    """
    Plot the phase plane for the deterministic model (Route I) 
    showing the vector field and the trajectory of the system.
    """
    X, Y = np.meshgrid(grid, grid)
    U, V = vector_field_det(model_type, grid, initial_conditions, params)
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


def gene_regulation_sde(y, t, dt,
                        a_a, a_b, b_a, b_b,
                        c_a, c_b, beta_a, beta_b,
                        gamma_a, gamma_b, n_a, n_b,
                        theta_a, theta_b, k_a, k_b,
                        delta_a, delta_b, m_a, m_b,
                        sigma_1a, sigma_2a, sigma_1b, sigma_2b,
                        noise, time_dependent_alpha):
    
    """
    Stochastic ODE model of gene regulation with Protein A inhibiting splicing of mRNA B,
    and Protein B promoting splicing of mRNA A (Route II). 
    This function includes noise terms for the transcription and splicing processes.
    """
    
    u_a, u_b, s_a, s_b, p_a, p_b = y

    # Noise terms 
    noise_terms = lambda sigma: sigma * np.random.normal(0, np.sqrt(dt)) if noise else 0

    # Time-dependent transcription rates 
    if time_dependent_alpha:
        alpha_a = c_a / (1 + np.exp(b_a * t - a_a))
        alpha_b = c_b / (1 + np.exp(b_b * t - a_b))
    else:
        # Use constant transcription rates (m_a, m_b)
        alpha_a = m_a
        alpha_b = m_b
    
    # Inhibition and promotion terms for splicing
    inhibition_splicing_b = 1 / (1 + (p_a / theta_a)**n_a)  # Inhibition of splicing of mRNA B by protein A
    promotion_splicing_a = (p_b / theta_b)**n_b / (1 + (p_b / theta_b)**n_b)  # Promotion of splicing of mRNA A by protein B
    
    # Transcription 
    du_a_dt = (alpha_a - beta_a * u_a) + noise_terms(sigma_1a)  # Transcription rate of pre-mRNA A
    du_b_dt = (alpha_b - beta_b * u_b) + noise_terms(sigma_1b)  # Transcription rate of pre-mRNA B
    
    # Splicing 
    ds_a_dt = (beta_a * u_a * promotion_splicing_a - gamma_a * s_a) + noise_terms(sigma_2a)  # Splicing of mRNA A
    ds_b_dt = (beta_b * u_b * inhibition_splicing_b - gamma_b * s_b) + noise_terms(sigma_2b)  # Splicing of mRNA B
    
    # Translation 
    dp_a_dt = k_a * s_a - delta_a * p_a  # Protein A production and degradation
    dp_b_dt = k_b * s_b - delta_b * p_b  # Protein B production and degradation
    
    return du_a_dt, du_b_dt, ds_a_dt, ds_b_dt, dp_a_dt, dp_b_dt


def solve_gene_regulation_sde(model_type, y0, t, params):
    """
    Solve the stochastic ODE model using Euler's method.
    This function simulates the time evolution of the system using a stochastic approach.
    """
    dt = t[1] - t[0]
    num_steps = len(t)
    sol = np.zeros((num_steps, len(y0)))
    sol[0] = y0

    for i in range(1, num_steps):
        dydt = model_type(sol[i-1], t[i-1], dt, **params)
        sol[i] = sol[i-1] + np.array(dydt) * dt

    return sol

def plot_mRNA_time_evolution_sde(t, solutions, save_path=None):
    """
    Plot the time evolution of mRNA concentrations for A and B (Route II).
    Error is shaded around the mean trajectory.
    """
    # Calculate mean and standard deviation across simulations
    solutions = np.array(solutions)
    mean_u_a, mean_u_b = np.mean(solutions[:, :, 0], axis=0), np.mean(solutions[:, :, 1], axis=0)
    mean_s_a, mean_s_b = np.mean(solutions[:, :, 2], axis=0), np.mean(solutions[:, :, 3], axis=0)
    std_u_a, std_u_b = np.std(solutions[:, :, 0], axis=0), np.std(solutions[:, :, 1], axis=0)
    std_s_a, std_s_b = np.std(solutions[:, :, 2], axis=0), np.std(solutions[:, :, 3], axis=0)

    fig, axs = plt.subplots(2, 1, figsize=(10, 12))

    # Subplot 1: u_a and u_b
    axs[0].plot(t, mean_u_a, label='Mean u_A', color='blue')
    axs[0].fill_between(t, mean_u_a - std_u_a, mean_u_a + std_u_a, color='blue', alpha=0.3, label='u_A ± std')
    axs[0].plot(t, mean_u_b, label='Mean u_B', color='red')
    axs[0].fill_between(t, mean_u_b - std_u_b, mean_u_b + std_u_b, color='red', alpha=0.3, label='u_B ± std')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_ylabel('Concentration')
    axs[0].set_title('Stochastic Gene Regulation: Transcription (u_A and u_B)')
    axs[0].legend()
    axs[0].grid()

    # Subplot 2: s_a and s_b
    axs[1].plot(t, mean_s_a, label='Mean s_A', color='green')
    axs[1].fill_between(t, mean_s_a - std_s_a, mean_s_a + std_s_a, color='green', alpha=0.3, label='s_A ± std')
    axs[1].plot(t, mean_s_b, label='Mean s_B', color='orange')
    axs[1].fill_between(t, mean_s_b - std_s_b, mean_s_b + std_s_b, color='orange', alpha=0.3, label='s_B ± std')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('Concentration')
    axs[1].set_title('Stochastic Gene Regulation: Splicing (s_A and s_B)')
    axs[1].legend()
    axs[1].grid()

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()

def vector_field_sde(model_type, grid, initial_conditions, params):
    """
    Calculate the vector field for the stochastic model (Route II).
    This function computes the derivatives of the protein concentrations
    at each point in the grid to plot the phase plot.
    """
    U = np.zeros((len(grid), len(grid)))
    V = np.zeros((len(grid), len(grid)))

    u_A, u_B, s_A, s_B = initial_conditions[0], initial_conditions[1], initial_conditions[2], initial_conditions[3]

    for i, p_A in enumerate(grid):
        for j, p_B in enumerate(grid):
            y0 = [u_A, u_B, s_A, s_B, p_A, p_B]
            dydt = model_type(y0, 0, 0.01, **params)
            U[j, i] = dydt[4]  # dp_a/dt
            V[j, i] = dydt[5]  # dp_b/dt

    return U, V


def plot_phase_plane_sde(model_type, initial_conditions, var_indices, grid, params, solution=None, save_path=None):
    """
    Plot the phase plane for the stochastic model (Route II)
    showing the vector field and the trajectory of the system.
    """
    X, Y = np.meshgrid(grid, grid)
    U, V = vector_field_sde(model_type, grid, initial_conditions, params)
    speed = np.sqrt(U**2 + V**2)  # Magnitude of the vector field

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
    Main function to run the gene regulation models and Viterbi algorithm.
    """
    # Viterbi setup
    states = ['Exon', 'Intron']
    observations = ['A', 'G', 'C', 'G', 'C']
    start_probs = {'Exon': 0.5, 'Intron': 0.5}
    trans_probs = {'Exon': {'Exon': 0.9, 'Intron': 0.1}, 'Intron': {'Exon': 0.2, 'Intron': 0.8}}
    emit_probs = {'Exon': {'A': 0.25, 'U': 0.25, 'G': 0.25, 'C': 0.25}, 'Intron': {'A': 0.4, 'U': 0.4, 'G': 0.05, 'C': 0.15}}

    prob, path = viterbi_algorithm(observations, states, start_probs, trans_probs, emit_probs)
    print("Most probable path:", path)
    print("Probability of the path:", prob)

    # Deterministic model
    t = np.linspace(0, 100, 500)
    initial_det = [0.8, 0.8, 0.8, 0.8]
    det_params = {
        'm_a': 2.35, 'm_b': 2.35,
        'n_a': 3, 'n_b': 3,
        'theta_a': 0.21, 'theta_b': 0.21,
        'gamma_a': 1.0, 'gamma_b': 1.0,
        'k_a': 1.0, 'k_b': 1.0,
        'delta_a': 1.0, 'delta_b': 1.0
    }

    sol_det = solve_gene_regulation_det(gene_regulation_det, initial_det, t, det_params)
    plot_mRNA_time_evolution_det(t, sol_det)

    grid = np.linspace(0, 4, 30)
    plot_phase_plane_det(
        model_type=gene_regulation_det,
        initial_conditions=initial_det,
        var_indices=[2, 3],  # indices for p_a and p_b
        grid=grid,
        params=det_params,
        solution=sol_det
    )

    # Stochastic model
    initial_sde = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8] # [u_a, u_b, s_a, s_b, p_a, p_b]

    stoch_params = {
    # Time-dependent transcription parameters
    'a_a': 1.0, 'a_b': 0.25,
    'b_a': 0.0005, 'b_b': 0.0005,
    'c_a': 2.0, 'c_b': 0.5,
    # Constant transcription rates (when not using time-dependent)
    'm_a': 2.35, 'm_b': 2.35,
    # Splicing and degradation rates
    'beta_a': 2.35, 'beta_b': 2.35,
    'gamma_a': 1.0, 'gamma_b': 1.0,
    # Hill function parameters
    'n_a': 3, 'n_b': 3,
    'theta_a': 0.21, 'theta_b': 0.21,
    # Translation and protein degradation
    'k_a': 1.0, 'k_b': 1.0,
    'delta_a': 1.0, 'delta_b': 1.0,
    # Noise parameters
    'sigma_1a': 0.05, 'sigma_2a': 0.05,
    'sigma_1b': 0.05, 'sigma_2b': 0.05,
    # Control flags
    'time_dependent_alpha': False,  # Set to True to use sigmoid function
    'noise': True  # Set to False for deterministic solution
}

    simulations = 100
    sde_results = [
        solve_gene_regulation_sde(gene_regulation_sde, initial_sde, t, stoch_params)
        for _ in range(simulations)
    ]

    plot_mRNA_time_evolution_sde(t, sde_results)
    plot_phase_plane_sde(
        model_type=gene_regulation_sde,
        initial_conditions=initial_sde,  # [u_a, u_b, s_a, s_b, p_a, p_b]
        var_indices=[4, 5],  # indices for p_a and p_b
        grid=grid,
        params=stoch_params,
        solution=sde_results[0]  # Mean trajectory across simulations
    )


if __name__ == "__main__":
    main()
