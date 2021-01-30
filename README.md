# StdErrorsCalibratedParams_Matlab
Code for "Standard Errors for Calibrated Parameters"

Some programs require the cvx package to be installed.
Installation instructions can be found here
http://cvxr.com/cvx/doc/install.html


Directories and Top-Level Programs
----------------------------------
Main/           Code for solving for the worst-case optimal variance
                matrix via semidefinite programming or quantile
                regression (as described in paper) and using these
                programs to test hypotheses about parameters in
                classical minimum distance settings.

Supporting/     Supporting code called by programs in Main/ to do
                standard procedures not directly related to the worst
                case variance problem (e.g. quantile regression, GMM,
                setting up parallel pools, finite differences, etc.).
                Anything in this folder is a candidate for swapping out
                for more efficient or standardized routines later on.

Test/           Code to test out the programs in Main/. The programs in
                this folder are all called by running Test_All.m from
                the top level directory of this repo.


Test_All.m      Run from main directory to execute all test examples in
                the Test/ directory, which includes a test that cvx is
                installed and running properly.


