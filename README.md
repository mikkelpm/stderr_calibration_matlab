# Standard Errors for Calibrated Parameters

Matlab function that computes worst-case standard errors (SE) for minimum distance estimators, given knowledge of only the marginal variances (but not correlations) of the matched moments.

The computed worst-case SE for the estimated parameters are sharp upper bounds on the true SE (which depend on the unknown moment correlation structure). For over-identified models, the package also computes the efficient moment selection that minimizes the worst-case SE. Additionally, the package can carry out tests of parameter restrictions or over-identifying restrictions.

**Reference:**
Cocci, Matthew D., and Mikkel Plagborg-Møller (2021), "Standard Errors for Calibrated Parameters", [arXiv:2109.08109](https://arxiv.org/abs/2109.08109)

Tested in: Matlab R2021a on Windows 10 PC (64-bit)

Other versions: [Python](https://github.com/mikkelpm/stderr_calibration_python)

## Contents

- [example.m](example.m): Simple example illustrating the main functionality of the package step by step

- [@MinDist](@MinDist): Matlab class for minimum distance estimation, standard errors, and testing

- [application/run_all.m](application/run_all.m): Application and simulation study based on [Alvarez & Lippi (2014)](http://dx.doi.org/10.3982/ECTA10662), see top of file for data requirements

- [tests](tests): Unit tests

## Requirements

For some functionality (such as joint testing) it is necessary to install the [cvx](http://cvxr.com/cvx/doc/install.html) Matlab package.

## Acknowledgements

This material is based upon work supported by the NSF under Grant #1851665. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the NSF.
