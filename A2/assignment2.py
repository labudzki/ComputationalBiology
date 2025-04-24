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
import os


# Verterbi algorithm for hidden Markov model
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


# Determistic model of gene regulation (Route I)
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
    plt.xlabel('Time (s)', fontsize=16)
    plt.ylabel('Concentration', fontsize=16)
    plt.title('Deterministic Gene Regulation', fontsize=20)
    plt.legend(fontsize=14)
    plt.grid()
    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()


def plot_phase_plane_det(model_type, initial_conditions, t, params, save_path=None):
    """
    Solve the full ODE system and plot the phase trajectory (p_a vs p_b).
    """
    solution = solve_gene_regulation_det(model_type, initial_conditions, t, params)
    p_a = solution[:, 2]
    p_b = solution[:, 3]
    
    plt.figure(figsize=(8, 6))
    plt.plot(p_a, p_b, lw=2, color='darkblue')
    plt.xlabel(r'Protein A ($p_a$)', fontsize=16)
    plt.ylabel(r'Protein B ($p_b$)', fontsize=16)
    plt.title('Protein Phase Trajectory (p_a vs p_b)', fontsize=20)
    plt.grid(True)

    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()


# Stochastic model of gene regulation (Route II)
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
    inhibition_splicing_b = (theta_a**n_a) / (theta_a**n_a + p_a**n_a) # Inhibition of splicing of mRNA B by protein A
    promotion_splicing_a = (p_b**n_b) / (theta_b**n_b + p_b**n_b)  # Promotion of splicing of mRNA A by protein B
    
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

def plot_mRNA_time_evolution_sde(t, num_simulations, initial_conditions, sde_params, save_path=None):
    """
    Plot the time evolution of mRNA concentrations for A and B (Route II).
    Individual simulations are plotted as semi-transparent lines.
    """
    # Run simulations and collect results
    solutions = []
    for _ in range(num_simulations):
        sol = solve_gene_regulation_sde(gene_regulation_sde, initial_conditions, t, sde_params)
        solutions.append(sol)

    solutions = np.array(solutions)

    # Extract individual trajectories
    u_a_vals, u_b_vals = solutions[:, :, 0], solutions[:, :, 1]
    s_a_vals, s_b_vals = solutions[:, :, 2], solutions[:, :, 3]

    fig, axs = plt.subplots(2, 1, figsize=(10, 12))

    # Plot transcription trajectories (u_A and u_B)
    for i in range(num_simulations):
        axs[0].plot(t, u_a_vals[i], color='blue', alpha=0.5, lw=0.8)
        axs[0].plot(t, u_b_vals[i], color='red', alpha=0.5, lw=0.8)
    axs[0].set_xlabel('Time (s)', fontsize=16)
    axs[0].set_ylabel('Concentration', fontsize=16)
    axs[0].set_title('Stochastic Gene Regulation: Transcription (u_A and u_B)', fontsize=20)
    axs[0].legend(['u_A', 'u_B'], fontsize=14)
    axs[0].grid()

    # Plot splicing trajectories (s_A and s_B)
    for i in range(num_simulations):
        axs[1].plot(t, s_a_vals[i], color='green', alpha=0.5, lw=0.8)
        axs[1].plot(t, s_b_vals[i], color='orange', alpha=0.5, lw=0.8)
    axs[1].set_xlabel('Time (s)', fontsize=16)
    axs[1].set_ylabel('Concentration', fontsize=16)
    axs[1].set_title('Stochastic Gene Regulation: Splicing (s_A and s_B)', fontsize=20)
    axs[1].legend(['s_A', 's_B'], fontsize=14)
    axs[1].grid()

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()


def plot_phase_plane_sde(sde_solutions, num_simulations, t, save_path):
    """
    Function to plot the phase plane for multiple stochastic gene regulation simulations.
    """
    sde_solutions = np.array(sde_solutions)  # Convert list of solutions into a numpy array
    p_a_vals = sde_solutions[:, :, 4]  # Extract protein A values from the solutions
    p_b_vals = sde_solutions[:, :, 5]  # Extract protein B values from the solutions

    # Plot the individual trajectories
    plt.figure(figsize=(8, 6))
    for i in range(num_simulations):
        plt.plot(p_a_vals[i], p_b_vals[i], color='grey', alpha=0.5)  # Light lines for individual trajectories

    # Plot the average trajectory in dark color
    avg_p_a = np.mean(p_a_vals, axis=0)
    avg_p_b = np.mean(p_b_vals, axis=0)
    plt.plot(avg_p_a, avg_p_b, color='darkblue', label='Average trajectory')

    plt.xlabel(r'Protein A ($p_a$)', fontsize=16)
    plt.ylabel(r'Protein B ($p_b$)', fontsize=16)
    plt.title('Stochastic Gene Regulation: Phase Trajectory (p_a vs p_b)', fontsize=20)
    plt.grid(True)
    plt.legend(fontsize=14)

    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()

def main():
    """
    Main function to run the gene regulation models and Viterbi algorithm.
    """
    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
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
    plot_mRNA_time_evolution_det(t, sol_det, save_path=os.path.join(results_dir, "route1_mRNA_dynamics.png"))

    plot_phase_plane_det(
        model_type=gene_regulation_det,
        initial_conditions=initial_det,
        t=t,
        params=det_params,
        save_path=os.path.join(results_dir, "route1_phase_plane.png")
    )

    # Stochastic model setup
    initial_sde = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
    sde_params = {
        'a_a': 1, 'a_b': 0.25,
        'b_a': 0.0005, 'b_b': 0.0005,
        'c_a': 2, 'c_b': 0.5,
        'beta_a': 2.35, 'beta_b': 2.35,
        'gamma_a': 1.0, 'gamma_b': 1.0,
        'n_a': 3, 'n_b': 3,
        'theta_a': 0.21, 'theta_b': 0.21,
        'k_a': 1.0, 'k_b': 1.0,
        'delta_a': 1.0, 'delta_b': 1.0,
        'm_a': 2.35, 'm_b': 2.35,
        'sigma_1a': 0.05, 'sigma_2a': 0.05,
        'sigma_1b': 0.05, 'sigma_2b': 0.05,
        'noise': True,
        'time_dependent_alpha': True
    }

    num_simulations = 20
    sde_solutions = []

    for _ in range(num_simulations):
        sol = solve_gene_regulation_sde(gene_regulation_sde, initial_sde, t, sde_params)
        sde_solutions.append(sol)

    plot_mRNA_time_evolution_sde(
        t,
        num_simulations,
        initial_sde,
        sde_params,
        save_path=os.path.join(results_dir, "route2_mRNA_dynamics.png")
    )

    plot_phase_plane_sde(
        sde_solutions,
        num_simulations,
        t,
        save_path=os.path.join(results_dir, "route2_phase_plane.png")
    )


if __name__ == "__main__":
    main()
