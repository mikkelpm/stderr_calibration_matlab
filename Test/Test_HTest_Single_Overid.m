function [] = Test_HTest_Single_Overid()

clearvars
rng(314)

%% Tuning params

  start_all_at_theta0 = true;
  fullyFeasible       = false;

%% Setup

  p    = 3; % Number of moments
  k    = 2; % Number of parameters
  Nrep = 200;

  % Set G matrix
  G = [1  0;
       1  1;
       1 -1];

  % Given G matrix, create h (moment) function
  h            = @(theta) G*theta;
  estimateFcn_ = @(theta, W, muhat) fminsearch(@(theta_) (h(theta_)-muhat)'*W*(h(theta_)-muhat), theta);
  GFcn_        = @(theta) ComputeG(h, theta);

  % Set true parameter values theta0
  theta0 = randn(k,1);

  % Compute true mu0
  mu0 = h(theta0);

  % Compute mu0 for overid test
  mu1 = mu0 + 2*ones(p,1);


  % Set nonsingular V matrix for empirical moments

    % True V matrix
    rho = 0.9;
    rho = 0;
    rho = 0.5;
    V   = [1 rho rho; rho 1 rho; rho rho 1];
    V   = diag([3 2 1])*V*diag([3 2 1]);

  % Generate muhats

    muhats = mvnrnd(mu0, V, Nrep)';


%% Pass through code to do testing

  I       = eye(k);

  % Loop over reps
  approaches  = {'naive', 'wcopt', 'wcopt_onestep', 'opt', 'opt_onestep'};
  Napproaches = length(approaches);
  for n = 1:Nrep
    fprintf('Simulation %d...\n', n)

    % Loop over parameters to test
    for l = 1:k
      res(n,l) = HTest_CompareEstimators(muhats(:,n), V, 1, I(:,l), theta0, true, theta0, h, @(theta,W) estimateFcn_(theta,W,muhats(:,n)), GFcn_, fullyFeasible);

      % Add J-test output
      for a = 1:Napproaches
        approach = approaches{a};
        fprintf('Adding overid test results for approach: %s...\n', approach);
        res(n,l).(approach) = HTest_Overid(res(n,l).(approach), 1, p, 0.05, strcmp(approach, 'wcopt'));
      end
    end
  end



%% Tables

  approaches  = fields(res);
  Napproaches = length(approaches);

  % Do coverage tables
  coverage_CI               = nan(Napproaches, k);
  rej_overid_joint          = nan(Napproaches, k);
  rej_overid_single         = nan(p, k, Napproaches);

  for a = 1:Napproaches
    fprintf('Rejection tables for approach %s...\n', approaches{a})
    for l = 1:k
      coverage_CI(a,l)         = nanmean(cellfun(@(s) s.coverage,          {res(:,l).(approaches{a})}));
      rej_overid_joint(a,l)    = nanmean(cellfun(@(s) s.overid_rej_joint,  {res(:,l).(approaches{a})}));
      rej_overid_single_       =         cellfun(@(s) s.overid_rej_single, {res(:,l).(approaches{a})}, 'un', 0);
      rej_overid_single(:,l,a) = nanmean(horzcat(rej_overid_single_{:}), 2);

    end
  end

  disp('Coverage of confidence intervals for params')
  disp('Expect to see at least 0.95 for everything except naive which uses wrong V in constructing std errs')
  coverage_CI

  disp('Joint overid test')
  disp('Expect to see at most 0.05 for everything except naive which uses wrong V in constructing std errs')
  disp('Rows indicate assumptions on V and W used to conduct test: naive, wcopt, wcopt_onestep, opt, opt_onestep')
  disp('Columns are which parameter we were estimating. Only relevant for wcopt case bc that changes W. Other approaches use same W, so should have same number across columsn')
  rej_overid_joint

  disp('Overid test of each moment individually')
  disp('Five separate matrices output')
  disp('Rows are the individual moments we matched')
  disp('Columns are which parameter we were estimating. Only relevant for wcopt case bc that changes W. Other approaches use same W, so should have same number across columsn')
  rej_overid_single






end
