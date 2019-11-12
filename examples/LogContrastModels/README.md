[![MATLAB FEX](https://img.shields.io/badge/MATLAB%20FEX-60536-green.svg)][fex]
[![Minimum Version](https://img.shields.io/badge/Requires-R2014a-blue.svg)][matlab]


Regression models for compositional data via log-contrast formulations
=========

Here, we consider the special but important case of estimating a log-contrast model for compositional covariate X
and continuous outcomes Y as shown below: 

<a href="https://www.codecogs.com/eqnedit.php?latex=Y=\log(X)\beta&space;&plus;&space;\sigma&space;\epsilon&space;\qquad&space;\text{s.t.}\qquad&space;C^T&space;\beta&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Y=\log(X)\beta&space;&plus;&space;\sigma&space;\epsilon&space;\qquad&space;\text{s.t.}\qquad&space;C^T&space;\beta&space;=&space;0" title="Y=\log(X)\beta + \sigma \epsilon \qquad \text{s.t.}\qquad C^T \beta = 0" /></a>

We consider joint estimation of regression vectors and scales using perspective M-estimation. The details about the optimization model, objective functions, and the proximal algorithms are found in [3]. 

The folder comprises code and data for reproducing the numerical experiments in [3]. 

  h = figure;




