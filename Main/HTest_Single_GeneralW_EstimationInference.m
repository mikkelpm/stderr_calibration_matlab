% Routine to estimate theta (given some W) and construct information
% necessary for inference about a single linear combination of
% parameters in theta under that W, within a CMD setting.
%
% W               Weighting matrix to use to compute thetahat
%
% V               V matrix for computing std errors. Must be filled in,
%                 no NaNs, since this function is going to compute
%                 standard errors given that V.
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
function [res] = HTest_Single_GeneralW_EstimationInference(W, V, h, lambda, Nobs, estimateFcn, GFcn)

  % Estimate parameter using given W
  thetahat = feval(estimateFcn, W);

  % Given parameter information, construct and package all other
  % information implied by that estimate which is necessary for
  % inference
  res = HTest_Single_GeneralW_InferenceInfo(thetahat, W, V, h, lambda, Nobs, GFcn);

end
