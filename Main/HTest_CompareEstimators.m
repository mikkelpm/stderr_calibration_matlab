% For Monte Carlo comparison of different testing approaches, when the
% truth is known.
%
% Inputs
% ------
% muhat                 Vector of estimated reduced-form moments
%
% V                     Full V matrix for computing std errors. Must be
%                       filled in.
%
% Nobs                  Number of obs, for properly normalizing stderrs
%
% lambda                Vector s.t. lambda*theta is the linear
%                       comb that we want to test a hypothesis about.
%
% theta0                Starting parameter value for optimization
%
% start_all_at_theta0   Whether to start all thetahat estimation
%                       routines at that theta0
%
% theta_true            True parameter value
%
% h                     Function of theta that gives model-implied moments
%
% estimateFcn_          Function that takes as input argument an initial
%                       theta (to start the optimization routine from)
%                       and a weighting matrix W, then runs an
%                       optimization routine to return the corresponding
%                       thetahat. Passed in this way so that the user
%                       can define the optimization routine outside of
%                       the function most efficiently given their
%                       problem.
%
% GFcn_                 Function that takes as input argument parameter
%                       value theta, then computes G(theta)
%                       corresponding to that theta. Passed in this way
%                       so that the user can define the numerical
%                       differentiation routine outside of the function
%                       most efficiently given their problem.
%
% fullyFeasible         If true, do fully feasible comparison of the
%                       estimators in the sense of using estimated G,
%                       i.e. G(thetahat) where thetahat is the estimated
%                       value. If false, use true G (i.e. G at
%                       theta_true).
%
function [res] = HTest_CompareEstimators(muhat, V, Nobs, lambda, theta0, start_all_at_theta0, theta_true, h, estimateFcn_, GFcn_, fullyFeasible)


%% General info and setup

  % Number of parameters and moments
  k = length(theta0);
  p = length(muhat);

  % Define some zero threshold, below which, something is "zero" up
  % to numerical precision. Used to identify "zero" eigenvalues, which
  % aren't numerically identical to zero.
  zero_thresh = 1e-8;

  % Structure to hold results
  res = struct();


%% Compute G at theta_true, which we may use in place of G at thetahat for testing purposes

  G.true = feval(GFcn_, theta_true);
  if rank(G.true) < size(G.true,2)
    disp('Warning: G of deficient rank')
    return
  end



%% Set the G fcn that will be used to compute the asymptotic variance in final expressions
% - When constructing asymptotic variance, need a consistent estor of G.
% - If fully feasible, will use final estimate of theta to construct G
% - Else use infeasible true G, no matter theta estimate.
% - This is for testing and sanity checks.

  if fullyFeasible
    GFcn = @(theta) feval(GFcn_, theta);
  else
    GFcn = @(theta) G.true;
  end


%% Estimate theta under naive W matrix, and get corresponding standard errors under

  % Naive
  % - Assumes (incorrectly) that all moments uncorrelated
  % - Then uses "efficient" weighting matrix under that assumption
  res.naive.name  = 'Naive';
  res.naive.V     = diag(diag(V));
  res.naive.W     = inv(res.naive.V);
  res.naive       = HTest_Single_GeneralW_EstimationInference(res.naive.W, res.naive.V, h, lambda, Nobs, @(W)feval(estimateFcn_, theta0, W), GFcn);



%% Whether we start subsequent optimizations from naive est or theta0
%% What to use as G, the truth or the naive estimate

  % Where to start subsequent optimizations from, the original theta0 or the
  % naive estimate
  if start_all_at_theta0
    theta0_use = theta0;
  else
    theta0_use = res.naive.theta;
  end
  estimateFcn = @(W)feval(estimateFcn_, theta0_use, W);

  % G to use in constructing asymptotic variance estimators
  % Note: Given the way GFcn derived above, if fullyFeasible=0, res.naive.G = G.true.
  Ghat0 = res.naive.G;

  % Construct fcn for getting onsetep estimator starting from naive est
  getOnestep_ = @(W) getOnestep(res.naive.h, W, Ghat0, res.naive.theta, muhat);


%% Compute wcopt W matrix, then est theta using it, then get standard errors


  % Worst-case optimal
  % - Assumes only diagonal
  % - Then estimates using worst-case optimal weighting matrix, both
  %   full optimization and one-step
  % - We have already done the first step in getting the naive
  %   estimates. Don't need to repeat. So just call Step2.
  res.wcopt.name                 = 'Worst-Case Optimal';
  res.wcopt_onestep.name         = 'Worst-Case Optimal, One-Step';
  [res.wcopt, res.wcopt_onestep] = HTest_Single_WorstCaseOptW_Step2(V, Ghat0, GFcn, lambda, Nobs, true, h, estimateFcn, true, muhat, res.naive.theta, res.naive.h);



%% Use full information optimal

  % Full-information optimal, one-step
  % - Requires the full V matrix is given/available
  % - Then W should be the optimal weighting matrix
  res.opt.name  = 'Full-Info Optimal';
  res.opt.V     = V;
  res.opt.W     = inv(V);
  res.opt       = HTest_Single_GeneralW_EstimationInference(res.opt.W, res.opt.V, h, lambda, Nobs, estimateFcn, GFcn);

  res.opt_onestep.name    = 'Full-Info Optimal, One-Step';
  res.opt_onestep.V       = res.opt.V;
  res.opt_onestep.W       = inv(res.opt.V);
  res.opt_onestep.theta   = getOnestep_(res.opt_onestep.W);
  res.opt_onestep         = HTest_Single_GeneralW_InferenceInfo(res.opt_onestep.theta, res.opt_onestep.W, res.opt_onestep.V, h, lambda, Nobs, GFcn);
  


%% Add muhat and coverage info

  approaches  = {'naive'; 'wcopt'; 'wcopt_onestep'; 'opt'; 'opt_onestep'};
  Napproaches = length(approaches);
  lcomb_true  = lambda'*theta_true;
  for n = 1:Napproaches
    approach = approaches{n};

    res.(approach).muhat      = muhat;
    res.(approach).lcomb_true = lcomb_true;
    res.(approach).coverage   = (lcomb_true > res.(approach).ci(1)) && (lcomb_true < res.(approach).ci(2));
    res.(approach).dist       = 100*abs(res.(approach).lcomb - lcomb_true)/lcomb_true;
  end



end




