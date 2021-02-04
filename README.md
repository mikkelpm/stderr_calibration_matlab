# StdErrorsCalibratedParams_Matlab
Code for "Standard Errors for Calibrated Parameters"

Some programs require the cvx package to be installed.
Installation instructions can be found here
http://cvxr.com/cvx/doc/install.html


TO-DO
-----
- Double check overid testing, where rejection rate for wcopt seems off
- Add nonlinear example to testing code


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




High Level Testing Functions Under Main/
----------------------------------------
HTest_CompareEStimators.m

  - For Monte Carlo testing of the different approaches.
  - Essentially calls everything below.
  - Requires known V


HTest_Single_GeneralW_EstimationInference.m

  - Get CI for single linear combination of params after a CMD
    estimation step. Steps:
    1. Estimate theta, given some W
    2. Construct CI for single linear combination of parameters in
       theta under that W, i.e. call HTest_Single_GeneralW_InferenceInfo.m

  - This is a self-contained estimation + inference routine (i.e. get
    the estimate for some W, build the info necessary for inference)
  - Does not implement two-step worst-case optimal aproach. Just one
    step of that.


HTest_Single_GeneralW_InferenceInfo.m

  - No estimation, just get inference info.
  - Compute valid CI (and other info) for single linear combination of
    parameters given estimates of theta (using some W).
  - Called by HTest_Single_GeneralW_EstimationInference.m and
    HTest_Single_WorstCaseOptW_Step2.m


HTest_Single_WorstCaseOptW_Full.m

  - Get CI from full two-step procedure for worst-case optimal
    inference about a single linear combination of parameters.
    Steps:

    1. Get initial consistent ests for theta and G given some W, i.e.
       call HTest_Single_GeneralW_EstimationInference.m

    2. Run second step that computes worst-case optimal W and
       reestimates, either by full reopt or 1-step iteration or both,
       i.e.  call HTest_Single_WorstCaseOptW_Step2.m


HTest_Single_WorstCaseOptW_Step2.m

  - Worst case computation + estimation + inference routine that returns
    a CI. Steps:

    1. Compute worst-case V and worst case optimal W given initial
       consistent estimate for G.

    2. Carry out a second step estimation of theta under that W, either
       via full reoptimization or one-step iteration.




Lower-Level Computations Under Main/
------------------------------------

ComputeWorstCaseOptW_Single.m

  - Compute worst-case optimal weighting matrix for testing hypotheses
    about a single parameter in classical minimum distance problem when
    only diag(V) known, using quantile regression approach outlined in
    paper.


ComputeCMDAvar.m

  - Formula to compute the asymptotic variance of a CMD estimator


getOnestep.m

  - Compute one-step estimator given some initial consistent estimates
    of theta and G and some W

SolveVarianceSDP

  - Solve for variance matrix in general SDP problem.
