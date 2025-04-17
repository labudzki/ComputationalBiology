#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 15:31:20 2025

@author: andrealabudzki
"""

import numpy as np

def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    path = {}

    for y in states:
        V[0][y] = start_p[y] * emit_p[y][obs[0]]
        path[y] = [y]

    for t in range(1, len(obs)):
        V.append({})
        newpath = {}

        for y in states:
            (prob, state) = max(
                [(V[t-1][y0] * trans_p[y0][y] * emit_p[y][obs[t]], y0) for y0 in states]
            )
            V[t][y] = prob
            newpath[y] = path[state] + [y]

        path = newpath

    (prob, state) = max([(V[-1][y], y) for y in states])
    return (prob, path[state])

states = ('Exon', 'Intron')
observations = ('A', 'G', 'C', 'G', 'C')
start_probability = {'Exon': 0.5, 'Intron': 0.5}
transition_probability = {
   'Exon' : {'Exon': 0.9, 'Intron': 0.1},
   'Intron' : {'Exon': 0.2, 'Intron': 0.8},
   }
emission_probability = {
   'Exon' : {'A': 0.25, 'U': 0.25, 'G': 0.25, 'C': 0.25},
   'Intron' : {'A': 0.4, 'U': 0.4, 'G': 0.05, 'C': 0.15},
   }

print(viterbi(observations, states, start_probability, transition_probability, emission_probability))
