# Attentional-Modes
Code associated with the manuscript "Dancing to the beats of the mammalian connectome: Asymmetric Dynamics Shape Rhythmic Brain Integration" (https://www.biorxiv.org/content/10.64898/2026.07.01.735524v1)

The study uses computationla models of the whole brain to examine how the structure of mammalian connectomes can generate periodic changes in information gating and alternating modes of attention.

refer to https://osf.io/24pmy for ./dependencies and ./figure_data

# Overview

This repository contains the code used to model slow coherence oscillations (SCOs) in connectome-based networks of Kuramoto oscillators, as described in the manuscript. Scripts expect their connectivity matrices, lookup tables, and precomputed phase timeseries in a "dependencies" directory and write outputs to a "figure_data" directory. Most contain human and macaque configurations that are toggled by hand.

"SCO_parameter_sweep.m" is the core simulation. It sweeps global coupling and conduction delay, recording metastability and the coherence timeseries of the whole cortex, temporo-parietal, and default mode networks. It takes a structural connectivity matrix and network lookup table as input, with human (Schaefer 1000) and macaque (regional map 82) paths provided. "AAL_validation.m" repeats this sweep on the 90-region AAL connectome.

"randomised_connectome.m" (edge-permuted connectomes) and "a2a_lorentz.m" (all-to-all network with Lorentzian frequencies) are the null models, and produce no periodic order parameter dynamics. "toy_hfreq.m" and "toy_hcoupling.m" are two-module toy models that sweep, respectively, intrinsic frequency and local coupling to validate the analytic beat expression.

"predicted_sco_frequencies.py" takes precomputed phase timeseries and, per network, reports the predicted SCO frequency (the difference between hemispheric synchronization frequencies) against the observed SCO frequency, writing a CSV table. "competitive_entrainment.m" takes a temporo-parietal phase timeseries and returns a wavelet transform of a sensory node alternating between network and stimulus entrainment. "gain_modulation.m" scales connectivity within one temporo-parietal hemisphere and measures the resulting SCO frequency across gain values, with and without delays.
