Regression models for compositional data via log-contrast formulations
=========

Here, we consider the special but important case of estimating a log-contrast model for compositional covariate X
and continuous outcomes Y that can also contain outliers o (aka the mean shift) as shown below: 

<a href="https://www.codecogs.com/eqnedit.php?latex=Y=\log(X)\beta&space;&plus;&space;o&space;&plus;&space;\sigma&space;\epsilon&space;\qquad&space;\text{s.t.}\qquad&space;C^T&space;\beta&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Y=\log(X)\beta&space;&plus;&space;o&space;&plus;&space;\sigma&space;\epsilon&space;\qquad&space;\text{s.t.}\qquad&space;C^T&space;\beta&space;=&space;0" title="Y=\log(X)\beta + o + \sigma \epsilon \qquad \text{s.t.}\qquad C^T \beta = 0" /></a>

We consider joint estimation of regression vectors and scales using perspective M-estimation. The details about the optimization model, objective functions, and the proximal algorithms are found in [[3]](https://arxiv.org/abs/1903.01050). 

The folder comprises code and data for reproducing the numerical experiments in [[3]](https://arxiv.org/abs/1903.01050). 

### Log-contrast models for BMI prediction from gut microbiome data ###


### Log-contrast models for pH prediction using soil microbiome data ###

The subsequent scripts describe the task of prediction soil pH from the relative abundances of 
soil microbiota. They reproduce the results from Section 4.2 in [[3]](https://arxiv.org/abs/1903.01050).

```Matlab
% Run the perspective log-contrast model with LS and Huber on pH data
% This script runs stability selection and stores the stability profiles
testLogContrastPHData
% 
% This script analyze stability selection results, does refitting, 
% and plots data used in Section 4.2 in [3]
analyzeStabSelPHData
```
This code reproduces the results in Appendix D in [[3]](https://arxiv.org/abs/1903.01050)
```Matlab
% Run the perspective log-contrast model with LS and Huber on pH data
% This script computes the entire regularization path 
testLogContrastPHTheo
% 
% This script analyze the theoretical lambda selection results,  
% and plots data used in Appendix D in [3]
analyzeLogContrastPHTheo
```

### Runtime analysis for BMI data ###

For the special case of constrained Lasso with joint scale estimation we compare the runtime of our
provably convergent scheme with a specialized coordinate descent scheme + line search scheme, 
proposed in [Shi et al., 2016](https://arxiv.org/abs/1603.00974). A modified version
of their code is available in the ```concomlasso.m ``` file.

```Matlab
% Run time test for the two algorithms on BMI data
runtimeBMI
```
This code reproduces Appendix C in [[3]](https://arxiv.org/abs/1903.01050).



