clear all;
addpath('../Main');
addpath('../Supporting');

% Tests of worst-case standard error functions


%% Settings

rng(202105053, 'twister'); % RNG seed
tol = 1e-6;                % Numerical tolerance


%% Analytical example in paper

for i=1:10 % Run several experiments...
    
    % Simulate random problem parameters
    G_elm = randn(4,1);
    G = [G_elm(1) 0; G_elm(2:3)'; 0 G_elm(4)];
    h = @(x) G*x;                   % Moment function
    theta = randn(size(G,2),1);     % True parameters
    mu = h(theta);                  % True moments
    sigma = exp(randn(size(mu)));   % Moment SE
    
    theta_estim_fct = @(W) (G'*W*G)\(G'*W*mu); % MD estimator given weight matrix W
    
    % Inefficient estimation/SE routine
    X0 = randn(length(mu));
    W0 = X0*X0'; % Ad hoc weight matrix
    res0 = WorstCaseSE(h, mu, sigma, theta_estim_fct, 'W', W0, 'eff', false);
    aux = (G'*W0*G)\(G'*W0);
    assert(norm(res0.r_theta_se-abs(aux)*sigma)<tol); % Check WCSE formula
    
    % Efficient estimation/SE routine
    res1 = WorstCaseSE(h, mu, sigma, theta_estim_fct, 'eff', true, 'one_step', true);
    res2 = WorstCaseSE(h, mu, sigma, theta_estim_fct, 'eff', true, 'one_step', false);
    
    % Check that we get same estimate from one-step as from full
    % optimization (due to linearity)
    assert(norm(res1.r_theta-res2.r_theta)<tol);
    
    % Estimate should also be the true parameter
    assert(norm(res1.r_theta-theta)<tol);
    
    % Check closed-form formula for solution (cf. Appendix)
    if sigma(1)*abs(G_elm(2)*G_elm(4)) <= sigma(2)*abs(G_elm(1)*G_elm(4))+sigma(3)*abs(G_elm(1)*G_elm(3))
        x_cf = [1/G_elm(1) 0 0]';
    else
        x_cf = [0 1/G_elm(2) -G_elm(3)/(G_elm(2)*G_elm(4))]';
    end
    assert(abs(res1.r_theta(1)-x_cf'*mu)<tol); % Estimate
    assert(abs(res1.r_theta_se(1)-abs(x_cf)'*sigma)<tol); % SE
    
end


%% cvx software

% Requires that software package cvx is installed
% http://cvxr.com/cvx/

assert(exist('cvx_begin', 'file')==2, 'cvx software is not installed - required for joint testing');

p = 5;

for i=1:10 % Run several experiments...

    % Simulate random parameters
    x = randn(p,1);
    G = x/(x'*x);
    mu = randn(p,1);
    sigma = exp(randn(p,1));
    
    % Worst-case SE
    res = WorstCaseSE(@(x) G*x, mu, sigma, [], 'theta', x'*mu, 'W', eye(p), 'eff', false);
    
    % Run cvx to find worst-case variance
    V = nan(p);
    V(eye(p)==1) = sigma.^2;
    [Vtilde, cvx_optval, cvx_status] = SolveVarianceSDP(x*x', V);
    V_tilde_corr = sqrt(diag(Vtilde)).\Vtilde./(sqrt(diag(Vtilde))');
    
    assert(strcmp(cvx_status, 'Solved')); % Check that cvx converged
    assert(norm(abs(V_tilde_corr(:))-1)<1e-4); % Check that worst-case V is perfect (+/-) correlation
    assert(abs(sqrt(cvx_optval)-res.r_theta_se)<tol); % Check that we get same WCSE

end
