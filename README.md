
PerspeCtive M-estimation (PCM) package 
=========

This is the PCM MATLAB package for perspective M-estimation. 
The package introduces an optimization model for maximum likelihood-type estimation (M-estimation) 
that generalizes a large class of known statistical models, including Huber’s concomitant M-estimation model, 
the scaled Lasso, \nu-Support Vector Machine Regression, and penalized estimation with structured sparsity. 
The model, termed perspective M-estimation, leverages the observation that a wide class of 
convex M-estimators with concomitant scale as well as structured norms are instances of perspective functions. 

The code builds on the following papers:

* [1] P. L. Combettes and C. L. Müller,[Perspective functions: Proximal calculus and applications in high-dimensional statistics](https://www.sciencedirect.com/science/article/pii/S0022247X16308071), J. Math. Anal. Appl., vol. 457, no. 2, pp. 1283–1306, 2018.
* [2] P. L. Combettes and C. L. Müller, [Perspective M-estimation via proximal decomposition](https://arxiv.org/abs/1805.06098), arXiv, 2018.
* [3] P. L. Combettes and C. L. Müller, [Regression models for compositional data: General log-contrast formulations, proximal optimization, and microbiome data applications](https://arxiv.org/abs/1903.01050), arXiv, 2019.

Developer: 
* Christian L. Müller, Simons Foundation (cmueller@flatironinstitute.org)

## Installation ##

The PCM package is self-contained. No external software needed. However, for testing the code base we rely 
on the [cvx package](http://cvxr.com/cvx/). 

After downloading the PCM package, use

```MATLAB
% This will add the folders to your MATLAB path
addPCM
```
to add all subfolders to your MATLAB path.

## Basic Usage ##

In the /examples/ folder you find several examples about the different modes of usage. 
Please refer to the README.md in the folder for further information.

## Extensions ##

