clear all;

% *** Standard errors for calibrated parameters: Example ***

% Consider a simple model with two structural parameters (theta_1,theta_2)
% and three reduced-form moments (mu_1,mu_2,mu_3).

% The theoretical mapping between parameters and moments is given by
% [mu_1; mu_2; mu_3] = [theta_1; theta_1+theta_2; 2*theta_2] = h(theta_1,theta_2)

% We observe the noisy estimates (hat{mu}_1,hat{mu}_2,hat{mu}_3) = (1.1,0.8,-0.1) of the true moments.

% The standard errors of the three empirical moments are (hat{sigma}_1,hat{sigma}_2,hat{sigma}_3)=(0.1,0.2,0.05).

% We will estimate the parameters (theta_1,theta_2) by minimum distance,
% matching the model-implied moments h(theta_1,theta_2) to the empirical moments:
% hat{theta} = argmin_{theta} (hat{mu}-h(theta))'*hat{W}*(hat{mu}-h(theta)).

% To compute standard errors for the estimated parameters, test hypotheses,
% and compute the efficient weight matrix hat{W}, we use the formulas in
% Cocci & Plagborg-Møller (2023) (https://arxiv.org/abs/2109.08109),
% which do not require knowledge of the correlation structure of the empirical moments.


%% Define the model

% We first import relevant packages and define the model and data.

% Define moment function h(.)
G = [1 0; 1 1; 0 2];
h = @(theta) G*theta;

% Define empirical moments and their s.e.
mu = [1.1,0.8,-0.1]';
sigma = [0.1,0.2,0.05]';

% Define MinDist object used in later analysis
obj = MinDist(h, mu, 'moment_se', sigma);

% (Note: In our simple example, we have a formula for the Jacobian of h(.) with respect to the parameters.
% This could be supplied to the "MinDist" call using the optional argument "moment_fct_deriv".
% The default behavior is to compute Jacobians numerically.)


%% Initial parameter estimates and standard errors

disp('*** INITIAL ESTIMATES ***');

% We first estimate the model using an ad hoc diagonal weight matrix
% hat{W}=diag(hat{sigma}_1^{-2},hat{sigma}_2^{-2},hat{sigma}_3^{-2}).
% The numerical optimization for computing the estimates (hat{theta}_1,hat{theta}_2) is started off at (0,0).

res = obj.fit('opt_init', zeros(2,1), 'eff', false); % eff=false: estimation based on ad hoc diagonal weight matrix

disp('Parameter estimates');
disp(res.estim');
disp('Standard errors');
disp(res.estim_se');

for i=1:2
    fprintf('%s%d\n', 'Worst-case moment var-cov matrix for estimating theta_', i);
    disp(res.worstcase_varcov{i});
end

% (Note 1: In this simple linear example, there exists a closed-form formula
% for the minimum distance estimator. This formula can be supplied to the "fit()"
% function using the optional argument "estim_fct".)

% (Note 2: In some cases the minimum distance parameter estimate may have 
% already been computed elsewhere. It can then be passed to the "fit()" function
% via the optional argument "param_estim". The function will compute the 
% corresponding standard errors without re-estimating the model.)


%% Test of parameter restrictions

disp('*** TEST OF PARAMETER RESTRICTIONS ***');

% Let us test whether the parameters theta_1 and theta_2 equal zero.

test_res = obj.test(res); % Tests are based on the "res" estimation results

disp('t-statistics for testing individual parameters');
disp(test_res.tstat');
disp('p-value of joint test');
disp(test_res.joint_pval);

% Using a 5% significance level, we cannot reject that theta_2 is zero
% individually based on its t-statistic. However, we can reject the joint 
% hypothesis that both parameters equal zero.

% Suppose we wanted to test the joint null hypothesis that (theta_1,theta_2)=(1,0).
% To do this, we first reformulate it as the hypothesis that the transformed vector
% r(theta_1,theta_2)=(theta_1-1,theta_2) has all elements equal to zero.
% We can then test the hypothesis as follows.

r = @(theta) theta-[1;0]; % Restriction function
res_restr = obj.fit('transf', r, 'opt_init', res.estim, 'eff', false); % Estimate the transformation r(theta)
test_res2 = obj.test(res_restr); % Test using the restriction function

disp('p-value of joint test');
disp(test_res2.joint_pval);


%% Over-identification test

disp('*** OVER-IDENTIFICATION TEST ***');

% Since we have more moments (3) than parameters (2), we can test the over-identifying restriction.
% One common way of doing this in applied work is to estimate the model 
% using only two of the moments and then checking whether the third, non-targeted moment
% at the estimated parameters is approximately consistent with the data.

weight_mat = diag([1/sigma(1)^2, 1/sigma(2)^2, 0]); % Weight matrix that puts no weight on third moment
res_justid = obj.fit('opt_init', zeros(2,1), 'eff', false, 'weight_mat', weight_mat);

disp('Just-identified parameter estimates');
disp(res_justid.estim');
disp('Model-implied moments');
disp(res_justid.moment_fit');

res_overid = obj.overid(res_justid); % Over-identification test based on just-identified estimates
disp('Error in matching non-targeted moment');
disp(res_overid.moment_error(3)); % The non-targeted moment is the third one
disp('Standard error');
disp(res_overid.moment_error_se(3));
disp('t-statistic');
disp(res_overid.tstat(3));

% Since the t-statistic is below 1.96, we can't reject the validity of the model at the 5% level.


%% Efficient estimation

disp('*** EFFICIENT ESTIMATES ***');

% The above estimation results relied on an ad hoc diagonal weight matrix.
% We can compute the weight matrix that minimizes the worst-case standard errors,
% and then report the corresponding estimates and standard errors.

res_eff = obj.fit('opt_init', zeros(2,1), 'eff', true); % Note: Efficient estimation (eff=true) is the default

disp('Efficient parameter estimates');
disp(res_eff.estim');
disp('Efficient standard errors');
disp(res_eff.estim_se');

for i=1:2
    fprintf('%s%d\n', 'Efficient moment loadings for estimating theta_', i);
    disp(res_eff.moment_loadings(:,i)');
end

% We see that theta_1 is estimated off the 1st moment only,
% while theta_2 is estimated off the 3rd moment only (up to small numerical error).

% (Note: The efficient estimates are not based on a single choice of weight matrix,
% since the efficient weight matrix depends on the specific parameter of interest.
% In the background, the analysis is actually run separately for each parameter.
% For this reason, it is not advised to use the "test()" or "overid()" commands
% with efficient estimation results. These commands are better used with estimation results
% that correspond to a single choice of weight matrix.)


%% Inference about transformed parameters

disp('*** TRANSFORMED PARAMETER ***');

% Suppose we want a confidence interval for the transformed parameter theta_1^2+theta_2.
% In a more realistic setting, this parameter might be some model-implied counterfactual of interest.
% We can do inference on transformed parameters using the "transf" argument
% to the "fit" function, as already used above.

res_transf = obj.fit('transf', @(theta) theta(1)^2+theta(2), 'opt_init', zeros(2,1)); % Efficient estimation (the default)

disp('Estimated transformation');
disp(res_transf.estim);
disp('Standard errors');
disp(res_transf.estim_se);

% (Note: If the gradient of the transformed parameter is available, 
% we can supply it to the "fit()" function using the optional "transf_deriv" argument.)


%% More information about the variance-covariance matrix

disp('*** MORE INFORMATION ABOUT VAR-COV MATRIX ***');

% Suppose we happen to also know that the first two empirical moments
% hat{mu}_1 and hat{mu}_2 are (asymptotically) independent.
% We can use this information to sharpen our inference about the parameters.
% First we define the known and unknown parts of the var-cov matrix of the empirical moments.

V = [sigma(1)^2 0 nan;
     0 sigma(2)^2 nan;
     nan nan sigma(3)^2];
disp('Var-cov matrix of moments');
disp(V);

% Then we define a "MinDist" object using this var-cov matrix and apply the estimation/testing routines.

obj_moreinfo = MinDist(h, mu, 'moment_varcov', V);
res_moreinfo = obj_moreinfo.fit('opt_init', zeros(2,1), 'eff', false);

disp('Initial estimates');
disp(res_moreinfo.estim');
disp('Standard errors');
disp(res_moreinfo.estim_se');

res_eff_moreinfo = obj_moreinfo.fit('opt_init', zeros(2,1), 'eff', true);
disp('Efficient estimates');
disp(res_eff_moreinfo.estim');
disp('Standard errors');
disp(res_eff_moreinfo.estim_se');


%% Full-information analysis

disp('*** FULL INFORMATION ***');

% Suppose finally that we know the entire var-cov matrix of the empirical moments. For example:

V_fullinfo = sigma .* [1, 0, 0.5; 0, 1, -0.7; 0.5, -0.7, 1] .* sigma';
disp('Var-cov matrix of moments');
disp(V_fullinfo);

% In this full-information setting, the econometric analysis is standard (Newey & McFadden, 1994).
% The estimation and testing routines work as before.

obj_fullinfo = MinDist(h, mu, 'moment_varcov', V_fullinfo);
res_fullinfo = obj_fullinfo.fit('opt_init', zeros(2,1), 'eff', false, 'weight_mat', diag(sigma.^(-2))); % Diagonal weight matrix

disp('Initial estimates');
disp(res_fullinfo.estim');
disp('Standard errors');
disp(res_fullinfo.estim_se');

res_eff_fullinfo = obj_fullinfo.fit('opt_init', zeros(2,1), 'eff', true); % Efficient weight matrix
disp('Efficient estimates');
disp(res_eff_fullinfo.estim');
disp('Standard errors');
disp(res_eff_fullinfo.estim_se');

test_res_fullinfo = obj_fullinfo.test(res_eff_fullinfo);
disp('p-value for joint test that both parameters are zero');
disp(test_res_fullinfo.joint_pval);

overid_res_fullinfo = obj_fullinfo.overid(res_eff_fullinfo);
disp('p-value for over-identification test');
disp(overid_res_fullinfo.joint_pval);
