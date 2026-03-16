# Spindle Synchro
A simple numerical integration routine and the corresponding analysis scripts for the model used in [Fiorenza et al., 2026](https://www.biorxiv.org/content/10.64898/2026.02.08.704347v1).

All code was developed and tested on MATLAB R2023b. No external packages or scripts are required for use. 

## Running simulations
This is done via `simulation.m`. Parameters are set to their nominal value by default. For a description of each, see Table S1. 

**Input:** none

**Output:** Equilibration time and steady-state values of various physical quantities, formatted as: 

`L | O_C | iKC | v_kf | v_bf (F_b)`

where `L` is spindle length, `O` is overlap length, and `iKC` is interkinetochore distance, all in microns. `v_kf` is k-fiber flux velocity and `v_bf` is bridging fiber flux velocity, both in microns per minute. Finally, `F_b` is the force experienced by the bridging fiber minus-end in piconewtons.

Additionally, two figures are produced by default. 
The first is a pseudo kymograph where the y-axis is time and the x-axis is position, and left/right centrosomes correspond to solid blue/orange lines, the plus-ends of left/right bridging fibers correspond to dashed green/light blue lines, and the plus-ends of left/right k-fibers correspond to dotted yellow/purple lines.
The second is a plot showing how all lengths, forces, and velocities change over time throughout the simulation. The force- and length-dependent scaling factors are also plotted for clarity, although they are technically redundant and could be obtained from the corresponding force and length plots.

To disable these plots, insert a block comment initializer `%{` at line 336 and 352 for the first and second figure, respectively.
### Parameter scans
To scan over 1 or 2 parameters, remove the block comment initializers `%{` at lines 93 and 549. By default, it is set up to explore `A_pk` and `A_m` parameters. 
Each loop iteration will output two text files: one with parameter values and one with steady-state values of all physical quantities. 
Create a folder which you would like these outputs to be stored in, then set `dir` equal to its location. 
It is recommended to disable the default plots (see above section) to reduce the time needed to complete the scan. 

## Recreating figure panels using analysis scripts
To create the synchronization plots which appear in Fig. 2, the `simulation.m` script itself was used. To recreate these plots, remove the block comments from lines 116 and 453. 
To disable the contribution from length-dependent effects, comment out the plus-end scaling factors in lines 182-185 or the minus-end scaling factors in lines 187-190. 

To create the contour and cross-sectional plots which appear in Figs. 3, 5, S1, S2, and S6, the `analysis_contours.m` script was used. 
To recreate these plots, perform a two-dimensional parameter scan of the parameters `A_pk` and `A_m`, `l_0__plus_k` and `L_0__minus`, or `v_0__plus_k_base` and `v_0__minus_base` and set `i_dir` in the analysis script to 2, 1, or 3, respectively. 
The script is set up for a 192 by 192 scan, which will be the output of the default for loop in `simulation.m` 

To create all other plots which appear in the manuscript, the `analysis_plots.m` script was used. Each block of code is labeled with the figure it was used to create. 
