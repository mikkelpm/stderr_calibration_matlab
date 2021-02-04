% Routine for testing a hypothesis about a single linear combination of
% parameters in CMD setting, with variance at least partially unknown.
%
% V               Var matrix with NaNs off diagonal reflectin unknown elements.
%                 TO DO: Update to use general SDP.
%
% Ghat0           What G to pass into the optimization problem.
%                 Since it may be costly to compute, this is passed in
%                 as an input argument rather that computed within this
%                 function.
%
% GFcn            Function that takes as input argument parameter value
%                 theta, then computes G(theta) corresponding to that
%                 theta. Passed in this way so that the user can define
%                 the numerical differentiation routine outside of the
%                 function most efficiently given their problem.
%
% lambda          lambda giving the linear combination lambda*theta that
%                 we want to test a hypothesis about.
%
% Nobs            Number of observations, for properly normalizing stderrs
%
% fullreopt       Whether to compute second step estimate of theta by
%                 full reoptimization. Both onestep and fullreopt can be
%                 true
%
% h               Function of theta that gives model-implied moments
%
% estimateFcn     Function that takes as input argument a weighing
%                 matrix W, then runs an optimization routine to return
%                 the corresponding thetahat. Passed in this way so that
%                 the user can define the optimization routine outside
%                 of the function most efficiently given their problem.
%
% onestep         Whether to compute second step estimate of theta by
%                 one-step iteration.
%
% muhat           Empirical moments, which are necessary for computing
%                 the one-step estimator
%
% thetahat        Initial consistent estimator to use in one-step update
%
% h_thetahat      h at that initial consistent estimator
%
function [res, res_onestep] = HTest_Single_WorstCaseOptW_Step2(V, Ghat0, GFcn, lambda, Nobs, fullreopt, h, estimateFcn, onestep, muhat, thetahat, h_thetahat)

  %% Compute worst-case optimal weighting matrix
  zero_thresh = 1e-8;
  [W_wcopt, stderr_check, V_wcopt, x_check, z_check] = ComputeWorstCaseOptW_Single(V, Ghat0, lambda, zero_thresh);


  %% Get final estimates for theta by fully reoptimizing
  if fullreopt
    % Do inference steps using this worst-case optimal weight matrix.
    % Then package some checks as well
    res = HTest_Single_GeneralW_EstimationInference(W_wcopt, V_wcopt, h, lambda, Nobs, estimateFcn, GFcn);
    res.stderr_check = stderr_check/sqrt(Nobs);
    res.x_check      = x_check;
    res.z_check      = z_check;
  else
    res = [];
  end


  %% Get final estimates by one-step iteration from initial estimates
  if ~isempty(onestep) && onestep
    res_onestep.theta       = getOnestep(h_thetahat, W_wcopt, Ghat0, thetahat, muhat);
    res_onestep             = HTest_Single_GeneralW_InferenceInfo(res_onestep.theta, W_wcopt, V_wcopt, h, lambda, Nobs, GFcn);
    res_onestep.lcomb_check = lambda'*thetahat - x_check'*(h_thetahat-muhat);

  else
    res_onestep = [];
  end

  return


end
