Regression models for compositional data via log-contrast formulations
=========

Here, we consider the special but important case of estimating a log-contrast model for compositional covariates X 
where each of the n rows comprises p-dimensional compositions (or relative abundances) and n continuous outcome variables 
Y that can also contain outliers o (in form of a (sparse) mean shift). The generative model thus reads: 

<a href="https://www.codecogs.com/eqnedit.php?latex=Y=\log(X)\beta&space;&plus;&space;o&space;&plus;&space;\sigma&space;\epsilon&space;\qquad&space;\text{s.t.}\qquad&space;C^T&space;\beta&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Y=\log(X)\beta&space;&plus;&space;o&space;&plus;&space;\sigma&space;\epsilon&space;\qquad&space;\text{s.t.}\qquad&space;C^T&space;\beta&space;=&space;0" title="Y=\log(X)\beta + o + \sigma \epsilon \qquad \text{s.t.}\qquad C^T \beta = 0" /></a>

We consider joint estimation of regression vectors &beta; and scales &sigma; using perspective M-estimation. The details about the optimization model, objective functions, and the proximal algorithms are found in [[3]](https://arxiv.org/abs/1903.01050). 

The folder comprises code and data for reproducing the numerical experiments in [[3]](https://arxiv.org/abs/1903.01050). 

### Log-contrast models for BMI prediction from gut microbiome data ###

The subsequent scripts describe the task of prediction body mass index (BMI) from the relative abundances of gut 
microbiota and two covariates from the COMBO dataset.

The first example applies log-contrast modeling on microbiome and covariate data with standard linear constraint.
```Matlab
% Run the perspective log-contrast model with LS and Huber on COMBO data
% This script runs stability selection, stores the stability profiles, and compute 
% the entire regularization path.
testLogContrastCOMBO
%
```
This analysis script reproduces Figure 3 in Section 4.1 in [[3]](https://arxiv.org/abs/1903.01050). 
```MATLAB
% This script analyze the regularization path selection results, does refitting, 
% and reproduces Figure 3 in Section 4.1 in [3].
analyzeLogContrastCOMBOTheo
```
This analysis script reproduces Figure 4 in Section 4.1 in [[3]](https://arxiv.org/abs/1903.01050). 
```MATLAB
% This script analyze the regularization path selection results, does refitting, 
% and reproduces Figure 4 in Section 4.1 in [3].
analyzeStabSelCOMBO
```
### Log-contrast models for BMI prediction from gut microbiome data with taxonomic coherence ###

The example applies log-contrast modeling on microbiome and covariate data with subcompostional constraint,
derived from phylogeny, i.e., OTUs are grouped by phylum and subcompositionally coherent with respect to this
grouping.

```MATLAB
% Run the perspective log-contrast model with LS and Huber on COMBO data.
% This script runs stability selection with and without subcompositional constraints 
% and stores the stability profiles
testLogContrastSubCompCOMBO
```
This analysis script reproduces (approximately) Figure 5 in Section 4.1 in [[3]](https://arxiv.org/abs/1903.01050). 

```MATLAB
% This script analyze the influence of the subcompositional constraint on stability selection
% and reproduces (up to fluctuations from random numbers) Figure 5 in Section 4.1 in [3].
analyzeStabSelDiffSubCompCOMBO

```

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



