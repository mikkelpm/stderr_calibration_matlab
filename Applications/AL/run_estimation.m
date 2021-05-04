clear all;
addpath('Functions');
addpath('../../Main');
addpath('../../Supporting');

% Application based on Alvarez & Lippi (ECMA 2014)


%% Settings

% Data
dat_file = 'wber.csv';  % CSV file with prices by store, week, and UPC (download from: https://www.chicagobooth.edu/research/kilts/datasets/dominicks)
select_store = 122;     % Single store ID to keep for the data analysis
drop_quantile = 0.99;   % Drop absolute price changes above this quantile

% Simulation
run_sim = true;             % Run simulation study?
numrep_sim = 1e4;           % Number of simulations
quant_grid_sim = linspace(sqrt(eps),1-sqrt(eps),1e3); % Grid for computing true quantiles of price changes (used for simulation)
alpha_sim = 0.05;           % Significance level of tests
rng_seed_sim = 202104151;   % Random number seed


%% Load and clean data

dat = readtable(dat_file);

% Keep only selected STORE
dat = dat(dat.STORE==select_store,:);

% Compute log price
dat = dat(~isnan(dat.PRICE),:); % Remove missing
dat = dat(dat.PRICE>=0.2 & dat.PRICE<=25, :); % Keep only items between 20 cents and 25 dollars

% Compute price changes
dat = sortrows(dat, {'UPC', 'WEEK'});
dat.lag_PRICE = [nan; dat.PRICE(1:end-1)]; % Lagged price
dat.lag_PRICE([true; diff(dat.UPC)~=0]) = nan; % Set lagged price to missing for first obs. of each UPC
dat.d_price = dat.PRICE-dat.lag_PRICE; % Absolute change
dat.dl_price = log(dat.PRICE)-log(dat.lag_PRICE); % Log change
dat.dl_price(abs(dat.d_price)<0.01) = 0; % Replace changes less than 1 cent by zero

% Drop extreme price changes (and NaN)
dl_price_quant = quantile(abs(dat.dl_price), drop_quantile);
dat = dat(abs(dat.dl_price)<=dl_price_quant,:);

% Non-zero price changes
dat.change = (dat.dl_price~=0);

% Center price changes (A&L model assumes symmetric distribution)
dat.dl_price_center = dat.dl_price-mean(dat.dl_price(dat.change));


%% Compute and display summary statistics

disp('SUMMARY STATISTICS');
disp(' ');

% Sample size
n = size(dat,1);
disp('Sample size');
disp(n);

% UPCs
upc_count = groupcounts(dat, 'UPC');

disp('Number of UPCs');
disp(size(upc_count,1));

disp('Avg. obs. per UPC');
disp(mean(upc_count.GroupCount));

% Moments
[mu_hat, V_hat] = estimate_moments(dat.change, dat.dl_price_center);
sigma_hat = sqrt(diag(V_hat));
sample_kurtosis = mu_hat(3)/mu_hat(2)^2;

disp('Fraction of nonzero price changes per week');
disp(mu_hat(1));

disp('Moments of price changes: std, kurt');
disp([sqrt(mu_hat(2)) sample_kurtosis]);

disp('Moments of absolute price changes: mean');
disp(mu_hat(4));

disp('Estimated correlation matrix of sample moments');
disp(sigma_hat.\V_hat./sigma_hat');


%% Estimation

disp(' ');
disp('ESTIMATION AND TESTING');

estim = estimate_model(mu_hat, V_hat);

disp(' ');
disp('Just-identified estimates: #prod, vol, sqrt(menu cost)');
disp(estim.justid.theta_hat');
disp('Full-info SE');
disp(estim.justid.fullinfo_se');
disp('SE under independence');
disp(estim.justid.indep_se');
disp('Worst-case SE');
disp(estim.justid.wc_se');
disp('Ratio of SE: full-info/worst-case');
disp(estim.justid.fullinfo_se'./estim.justid.wc_se');
disp('Ratio of SE: independence/worst-case');
disp(estim.justid.indep_se'./estim.justid.wc_se');

disp(' ');
disp('Over-ID test moment errors: avg#changes, std, 4th, avg-abs');
disp(estim.overid_test.errors');
disp('Full-info SE');
disp(estim.overid_test.fullinfo_se');
disp('SE under independence');
disp(estim.overid_test.indep_se');
disp('Worst-case SE');
disp(estim.overid_test.wc_se');
disp('Ratio of SE: full-info/worst-case');
disp(estim.overid_test.fullinfo_se'./estim.overid_test.wc_se');
disp('Ratio of SE: independence/worst-case');
disp(estim.overid_test.indep_se'./estim.overid_test.wc_se');
disp('Full-info p-values');
disp(2*normcdf(-abs(estim.overid_test.errors./estim.overid_test.fullinfo_se)'));
disp('p-values under independence');
disp(2*normcdf(-abs(estim.overid_test.errors./estim.overid_test.indep_se)'));
disp('Worst-case p-values');
disp(2*normcdf(-abs(estim.overid_test.errors./estim.overid_test.wc_se)'));

disp(' ');
disp('Worst-case efficient estimates: #prod, vol, sqrt(menu cost)');
disp(estim.wceff.theta_hat');
disp('Worst-case SE');
disp(estim.wceff.se');
disp('Ratio of worst-case SE: efficient/just-ID');
disp(estim.wceff.se'./estim.justid.wc_se');

disp(' ');
disp('Full-information efficient estimates: #prod, vol, sqrt(menu cost)');
disp(estim.fullinfo.theta_hat');
disp('SE');
disp(estim.fullinfo.se');
disp('Ratio of SE: full-info/worst-case-efficient');
disp(estim.fullinfo.se'./estim.wceff.se');


%% Simulation study, using empirically estimated parameters

if ~run_sim
    return;
end

% Ordering of procedures:
% 1: just-ID, full-info SE
% 2: just-ID, SE under independence
% 3: just-ID, worst-case SE
% 4: all moments, efficient, full-info SE
% 5: all moments, worst-case efficient, worst-case SE

% Preliminaries
theta_sim = estim.justid.theta_hat; % True parameters in simulations
n_sim = n; % Sample size for simulations
[mu_sim, y_bar_sim] = moment_function(theta_sim); % True moments
quantiles_sim = quantile_price(theta_sim(1), y_bar_sim, quant_grid_sim); % Compute quantiles of price change distribution on grid (for quick simulation below)

% Run simulations
sim_theta_hat = nan(numrep_sim,3,5);
sim_se = sim_theta_hat;
sim_overid_tstat = nan(numrep_sim,3); % Last dimension: 1=efficient, 2=independent, 3=worst-case

rng(rng_seed_sim, 'twister');
sim_seeds = randi(2^32-1,numrep_sim,1); % RNG seeds for parallel workers

disp('Simulating...');
poolobj = parpool; % Start parallel pool
timer = tic;

parfor i=1:numrep_sim
% for i=1:numrep_sim
    
    % Simulate data
    rng(sim_seeds(i), 'twister');
    [change_sim, price_sim] = simulate_data(mu_sim, quant_grid_sim, quantiles_sim, n_sim);
    
    % Estimate moments
    [mu_hat_sim, V_hat_sim] = estimate_moments(change_sim, price_sim);
    
    % Estimate and test model
    try
        estim_sim = estimate_model(mu_hat_sim, V_hat_sim);
    catch ME
        fprintf('%s%d\n', 'Error for i=', i);
        warning(ME.msgtext);
        continue;
    end
    
    % Store estimates and SE
    sim_theta_hat(i,:,:) = [repmat(estim_sim.justid.theta_hat,1,3) estim_sim.fullinfo.theta_hat estim_sim.wceff.theta_hat];
    sim_se(i,:,:) = [estim_sim.justid.fullinfo_se estim_sim.justid.indep_se estim_sim.justid.wc_se estim_sim.fullinfo.se estim_sim.wceff.se];
    
    % Compute over-ID t-statistics
    sim_overid_tstat(i,:) = abs(estim_sim.overid_test.errors(4))./[estim_sim.overid_test.fullinfo_se(4) estim_sim.overid_test.indep_se(4) estim_sim.overid_test.wc_se(4)];
    
    % Print progress
    if mod(i,ceil(numrep_sim/50))==0
        fprintf('%s%3d%s\n', repmat(' ',1,round(50*i/numrep_sim)), round(100*i/numrep_sim), '%');
    end
    
end

sim_elapsed_time = toc(timer);
delete(poolobj);

disp('Done. Elapsed time (minutes):');
disp(sim_elapsed_time/60);

% Collect results
sim = struct;
sim.theta_hat = sim_theta_hat;
sim.se = sim_se;
sim.overid_tstat = sim_overid_tstat;
sim.seeds = sim_seeds;
sim.elapsed_time = sim_elapsed_time;
clearvars sim_theta_hat sim_se sim_overid_tstat sim_seeds sim_elapsed_time;

% Coverage/rejection rates
cv_sim = norminv(1-alpha_sim/2); % Normal critical value
sim.theta_hat_error = sim.theta_hat-permute(theta_sim, [3 1 2]); % Estimation error
sim.cover = (abs(sim.theta_hat_error./sim.se)<cv_sim); % CI coverage indicator
sim.length = 2*sim.se*cv_sim; % CI length
sim.overid_reject = sim.overid_tstat>cv_sim; % Over-ID rejection indicator

sim.bias = permute(mean(sim.theta_hat_error,1,'omitnan'), [2 3 1]);
sim.var = permute(var(sim.theta_hat,0,1,'omitnan'), [2 3 1]);
sim.mse = sim.bias.^2+sim.var;
sim.cover_rate = permute(mean(sim.cover,1,'omitnan'), [2 3 1]);
sim.avg_length = permute(mean(sim.length,1,'omitnan'), [2 3 1]);
sim.med_length = permute(median(sim.length,1,'omitnan'), [2 3 1]);
sim.overid_reject_rate = mean(sim.overid_reject,1,'omitnan');

% Display results
disp('RMSE relative to truth: #prod, vol, sqrt(menu cost)');
disp(sqrt(sim.mse')./theta_sim');

disp('CI coverage rate: #prod, vol, sqrt(menu cost)');
disp(sim.cover_rate');

disp('CI avg. length: #prod, vol, sqrt(menu cost)');
disp(sim.avg_length');

disp('Over-ID rejection rate: avg-abs');
disp(sim.overid_reject_rate');