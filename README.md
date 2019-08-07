
PerspeCtive M-estimation (PCM) package 
=========

PerspeCtive M-estimation (PCM) 

This is the PCM MATLAB package for perspective M-estimation accompanying the paper
[Proximal Analysis for Perspective M-estimation](...). The package introduces an optimization model 
for maximum likelihood-type estimation (M-estimation) that generalizes 
a large class of known statistical models, including Huber’s concomitant M- estimation model, 
the scaled Lasso, Support Vector Machine Regression, and penalized estimation with structured sparsity. 
The model, termed perspective M-estimation, leverages the observation that convex M-estimators with 
concomitant scale as well as structured norms are instances of perspective functions. 

The code developed here also builds on prior work:
[Perspective functions: Proximal calculus and applications in high-dimensional statistics](https://www.sciencedirect.com/science/article/pii/S0022247X16308071)

Authors: Patrick L. Combettes, North Carolina State University (plc@math.ncsu.edu),
Christian L. Mueller, Simons Foundation (cmueller@flatironinstitute.org)

Developer: Christian L. Mueller, Simons Foundation (cmueller@flatironinstitute.org)

## Installation ##

The package is mostly self-contained. No external software needed. Run the addPCM.m script to include 
the subfolders to your MATLAB path. Some of the functions for testing require an installation of cvx.

## Basic Usage ##

In the /examples/ folder you find several examples about the different modes of usage. 
Please refer to the README.md in the folder for further information.

## Extensions ##

