#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 15:31:20 2025

@author: andrealabudzki
"""

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


if __name__ == "__main__":
    main()
