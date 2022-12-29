# A simple, efficient, mass conservative approach to solving Richards’ Equation

Ireson, A.M.
andrew.ireson@usask.ca

This repo contains various models, mostly written in python, to reproduce the figures in Ireson et al. (2022, submitted), A simple, efficient, mass conservative approach to solving Richards’ Equation, GMD, submitted.

There are four scenarios modeled in this repo, all of which are described in the paper:

1. Celia's problem: generates Figures 2 and 3

2. Miller's problem: generates Figure 4

3. Mathias' problem: generates Figure 5

4. Infiltration problem: generates Figures 6 and 7

All models can be run using the Makefiles in the subfolders. In each case, it is necessary to first run each model, within the subfolders, and then in the folder "plot_result" there is a Makefile that will produce the figures from the paper

We also include a MATLAB implementation of openRE that solves the Infiltration problem in the folder `infiltrationproblem/modelbenchmarks/matlab`
