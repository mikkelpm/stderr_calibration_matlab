clear all;
addpath('../');

% Unit tests of minimum distance inference functions


% Define small-scale problem used for some tests

G = [1 0; 1 1; 0 2];
h = @(x) G*x;
theta = [1; 1];
mu = h(theta);
sigma = [1; 2; 0.5];

V_fullinfo = sigma.*[1,0.5,0.5; 0.5,1,0.5; 0.5,0.5,1].*sigma';

V_blockdiag = V_fullinfo;
V(1,2:3) = nan;
V(2:3,1) = nan;

tol = 1e-6; % Numerical tolerance


%% Closed-form formulas

obj = MinDist(h, mu, 'moment_se', sigma);

% Estimation with default weight matrix;
res = obj.fit('opt_init', zeros(size(theta)), 'eff', false);
W = diag(1./sigma.^2);
aux = W*G/(G'*W*G);
assert(norm(res.moment_loadings - aux)<tol);
assert(norm(res.estim_se - abs(aux)'*sigma)<tol);

% Efficient estimation (see formulas in paper appendix)
res_eff = obj.fit('opt_init', zeros(size(theta)));
if sigma(1)*abs(G(2,1)*G(3,2)) <= sigma(2)*abs(G(1,1)*G(3,2))+sigma(3)*abs(G(1,1)*G(2,2))
    x = [1/G(1,1); 0 ; 0];
else
    x = [0; 1/G(2,1); -G(2,2)/(G(2,1)*G(3,2))];
end
assert(norm(x'*mu - res_eff.estim(1))<tol);
assert(norm(abs(x)'*sigma - res_eff.estim_se(1))<tol);


%% Known diagonal

% Limited information
obj = MinDist(h, mu, 'moment_se', sigma);
[res, res_eff] = run_tests(obj, G, mu, theta, tol);
for i=1:length(theta)
    assert(norm(diag(res.worstcase_varcov{i}) - sigma.^2)<tol);
end
assert(all(res_eff.estim_se<=res.estim_se));

% Test that full optimization gives the same as one-step (due to linear moment function)
res2_eff = obj.fit('opt_init', zeros(size(theta)), 'one_step', false);
assert(norm(res_eff.estim - res2_eff.estim)<tol);
assert(norm(res_eff.estim_se - res2_eff.estim_se)<tol);

% Full information
obj_fullinfo = MinDist(h, mu, 'moment_varcov', V_fullinfo);
[res_fullinfo, res_eff_fullinfo] = run_tests(obj_fullinfo, G, mu, theta, tol);
assert(all(res_fullinfo.estim_se<=res.estim_se));
assert(all(res_eff_fullinfo.estim_se<=res_eff.estim_se));


%% Positive semidefinite problem

obj = MinDist(h, mu, 'moment_se', sigma);
res = obj.fit('opt_init', zeros(size(theta)), 'eff', false);
obj.diag_only = false; % Force PSD programming
res2 = obj.fit('opt_init', zeros(size(theta)), 'eff', false);
assert(norm(res2.estim_se - res.estim_se)<tol);
for i=1:length(theta)
    assert(norm(res2.worstcase_varcov{i} - res.worstcase_varcov{i})<tol);
end


%% Block diagonal

obj = MinDist(h, mu, 'moment_varcov', V_blockdiag);
res = obj.fit('opt_init', zeros(size(theta)), 'eff', false);
res_eff = obj.fit('opt_init', zeros(size(theta)));
obj.blockdiag_only = false; % Force PSD programming
res2 = obj.fit('opt_init', zeros(size(theta)), 'eff', false);
res2_eff = obj.fit('opt_init', zeros(size(theta)));

assert(all(res_eff.estim_se<=res.estim_se));
assert(norm(res.estim_se - res2.estim_se)<tol);
assert(norm(res_eff.estim_se - res2_eff.estim_se)<tol);


%% High-dimensional example

% Generate random high-dimensional problem with known block diagonal
blocks_num = 4; % Number of blocks; 1st block is 1x1, 2nd block is 2x2, etc.
k = 4; % Number of parameters
rng(123, 'twister'); % Random seed
G = randn(blocks_num*(blocks_num+1)/2,k);
h = @(x) G*x;
mu = G*(0:k-1)';
V_blocks = cell(1,blocks_num);
for i=1:blocks_num
    x = randn(i);
    V_blocks{i} = x*x';
end
V = blkdiag(V_blocks{:});
V(V==0) = nan;

% Estimate
obj = MinDist(h, mu, 'moment_varcov', V);
res = obj.fit('opt_init', zeros(k,1), 'eff', false);
res_eff = obj.fit('opt_init', zeros(k,1));
obj.blockdiag_only = false; % Force PSD programming
res2 = obj.fit('opt_init', zeros(k,1), 'eff', false);
res2_eff = obj.fit('opt_init', zeros(k,1));

assert(all(res_eff.estim_se<=res.estim_se));
assert(norm(res.estim_se - res2.estim_se)<tol);
assert(norm(res_eff.estim_se - res2_eff.estim_se)<0.002); % Allow higher tolerance
