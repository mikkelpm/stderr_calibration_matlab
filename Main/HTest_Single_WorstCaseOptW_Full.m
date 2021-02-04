% Implements full two-step procedure to worst-case optimal inference,
% returning a valid confidence interval.
%
% W0              Initial weighting matrix to use to compute thetahat in
%                 first step.
%
% V               V matrix for computing std errors, with NaNs on
%                 off-diagonal.
%                 TO DO: Allow non-diagonal V
%
% h               The structural function mapping parameters into moments
%
% lambda          Vector s.t. lambda*theta is the linear combination
%                 that we want to test a hypothesis about.
%
% Nobs            Number of observations, for properly normalizing stderrs
%
% estimateFcn     Function that takes as input argument a weighing
%                 matrix W, then runs an optimization routine to return
%                 the corresponding thetahat. Passed in this way so that
%                 the user can define the optimization routine outside
%                 of the function most efficiently given their problem.
%
% GFcn            Function that takes as input argument parameter value
%                 theta, then computes G(theta) corresponding to that
%                 theta. Passed in this way so that the user can define
%                 the numerical differentiation routine outside of the
%                 function most efficiently given their problem.
%
% fullreopt       Whether to compute second step estimate of theta by
%                 full reoptimization. Both onestep and fullreopt can be
%                 true
%
% onestep         Whether to compute second step estimate of theta by
%                 one-step iteration.
%
function [res, res_onestep] = HTest_Single_WorstCaseOpt_Full(W0, V, h, lambda, Nobs, estimateFcn, GFcn, fullreopt, onestep)

  % Step 1: Get initial consistent estimates under given W0
  res = HTest_Single_GeneralW_EstimationInference(W0, V, h, lambda, Nobs, estimateFcn, GFcn);

  % Step 2: Get final estimates using worst-case optimal W
  [res, res_onestep] = HTest_Single_WorstCase_Step2(V, res.G, GFcn, lambda, Nobs, fullreopt, h, estimateFcn, onestep, muhat, thetahat, h_thetahat);

end
