% Routine run after estimation of theta (given some W) to compute all
% information that's necessary for inference about a single linear
% combination of parameters in theta under that W.
%
% thetahat        Parameter estimate obtained from using W
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
% GFcn            Function that takes as input argument parameter value
%                 theta, then computes G(theta) corresponding to that
%                 theta. Passed in this way so that the user can define
%                 the numerical differentiation routine outside of the
%                 function most efficiently given their problem.
%
function [res] = HTest_Single_GeneralW_InferenceInfo(thetahat, W, V, h, lambda, Nobs, GFcn)

  % Implied moments under that theta
  h = feval(h, thetahat);

  % G matrix under that theta
  G = feval(GFcn, thetahat);


  % Get the target linear combination of parameters, which we are
  % testing a hypothesis about
  lcomb = lambda'*thetahat;

  % Vector x as defined in paper
  x = ( lambda'*((G'*W*G) \ (G'*W)) )';

  % Standard error for the linear combination
  stderr = sqrt(ComputeCMDAvar(G, W, V, lambda')/Nobs);

  % Confidence interval
  cv = 1.96;
  ci = lcomb +  [-cv cv]*stderr;

  % Package results
  res.theta  = thetahat;
  res.W      = W;
  res.V      = V;
  res.h      = h;
  res.G      = G;
  res.lcomb  = lcomb;
  res.x      = x;
  res.stderr = stderr;
  res.ci     = ci;

end
