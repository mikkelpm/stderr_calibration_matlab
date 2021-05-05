# Standard Errors for Calibrated Parameters

Matlab function that computes worst-case standard errors (SE) for minimum distance estimators, given knowledge of only the marginal variances (but not correlations) of the matched moments. The computed worst-case SE for the estimated parameters are sharp upper bounds on the true SE (which depend on the unknown moment correlation structure). For over-identified models, the function also computes the efficient moment selection that minimizes the worst-case SE.

**Reference:**
Cocci, Matthew D., and Mikkel Plagborg-Møller (2019), "Standard Errors for Calibrated Parameters", https://scholar.princeton.edu/mikkelpm/calibration

**Requirements:**
To perform joint hypothesis tests, it is necessary to install the [cvx](http://cvxr.com/cvx/doc/install.html) Matlab package. This package is not required for other functionality.

Tested in: Matlab R2021a on Windows 10 PC (64-bit)

## Contents

**[Main](Main):** main Matlab functions
- [WorstCaseSE.m](Main/WorstCaseSE.m): main function for computing worst-case SE and efficient estimates (can also do full-information estimation)

**[Applications](Applications):** empirical examples
- [AL/run_estimation.m](Applications/AL/run_estimation.m): application and simulation study based on Alvarez & Lippi (2014), see file for data requirements

**[Supporting](Supporting):** auxiliary Matlab functions

**[Test](Test):** unit tests

## Example

``` Matlab
addpath('Main', 'Supporting');

% Minimum distance problem (simple example with p=3, k=2)
h = @(theta) [1 0; 1 1; 0 2]*theta; % Moment function h(theta)
mu = [1; 2; 2];                     % Empirically estimated moments
sigma = [1; 2; 0.5];                % SE of empirical moments
theta_estim_fct = @(W) fminunc(@(theta) (mu-h(theta))'*W*(mu-h(theta)), zeros(2,1));
  % Estimator as a function of weight matrix
  % (closed form is available in this case, but we ignore that here)

% Compute efficient worst-case estimate and SE
res = WorstCaseSE(h, mu, sigma, theta_estim_fct);
disp('theta estimate');
disp(res.r_theta');
disp('SE');
disp(res.r_theta_se');
```
See the top of the file [WorstCaseSE.m](Main/WorstCaseSE.m) for additional optional inputs and a detailed description of the output structure. See the [empirical application](Applications/AL/run_estimation.m) for further illustrations of the functionality (e.g., over-identification tests and parameter restriction tests).


