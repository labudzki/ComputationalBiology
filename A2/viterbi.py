#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 15:31:20 2025

@author: andrealabudzki
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def viterbi_algorithm(obs, states, start_probs, trans_probs, emit_probs):
    """
    Viterbi algorithm for finding the most probable sequence of hidden states
    given a sequence of nucleotides.
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

def gene_regulation_model_1(y, t, m_a, m_b, p_a, p_b, n_a, n_b, theta_a, theta_b, gamma_a, gamma_b, r_a, r_b, k_a, k_b, delta_a, delta_b): 
    """
    Gene regulation model for mRNA and protein concentrations for route I. 
    """
    dr_a_dt = m_a * (p_b**n_b)/(p_b**n_b + theta_b**n_b) - gamma_a * r_a
    dr_b_dt = m_b * (theta_a**n_a)/(p_a**n_a + theta_a**n_a) - gamma_b * r_b

    dp_a_dt = k_a * r_a - delta_a * p_a
    dp_b_dt = k_b * r_b - delta_b * p_b
    return dr_a_dt, dr_b_dt, dp_a_dt, dp_b_dt

def gene_regulation_model_2(y, t, m_a, m_b, p_a, p_b, n_a, n_b, theta_a, theta_b, gamma_a, gamma_b, r_a, r_b, k_a, k_b, delta_a, delta_b):
    """
    placeholder for the second gene regulation model (route II)
    """

def solve_gene_regulation_model(gene_regulation_model, initial_conditions, t, params):
    """
    Solves the specified gene regulation model using the odeint solver.
    """
    solution = odeint(gene_regulation_model, initial_conditions, t, args=(
        params['m_a'], params['m_b'], params['p_a'], params['p_b'],
        params['n_a'], params['n_b'], params['theta_a'], params['theta_b'],
        params['gamma_a'], params['gamma_b'], params['r_a'], params['r_b'],
        params['k_a'], params['k_b'], params['delta_a'], params['delta_b']
    ))
    return solution

def plot_mRNA_time_evolution(t, solution):
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
    plt.show()

def main():
    """
    Main function to run the Viterbi algorithm with given parameters.
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
    t = np.linspace(0, 10000, 1000)  

    params = {
        'm_a': 2.35, 'm_b': 2.35,  # Max transcription rates
        'p_a': 0.8, 'p_b': 0.8,    # Initial protein concentrations
        'n_a': 3, 'n_b': 3,        # Hill coefficients
        'theta_a': 0.21, 'theta_b': 0.21,  # Binding thresholds
        'gamma_a': 1.0, 'gamma_b': 1.0,    # mRNA degradation rates
        'r_a': 0.8, 'r_b': 0.8,    # Initial mRNA concentrations
        'k_a': 1.0, 'k_b': 1.0,    # Translation rates
        'delta_a': 1.0, 'delta_b': 1.0  # Protein degradation rates
    }

    solution = solve_gene_regulation_model(
        gene_regulation_model_1,
        initial_conditions,
        t,
        params
    )

    plot_mRNA_time_evolution(t, solution)



if __name__ == "__main__":
    main()
