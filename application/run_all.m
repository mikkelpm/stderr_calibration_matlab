clear all;
addpath('functions');
addpath('../');

% Application based on Alvarez & Lippi (ECMA 2014)


%% Settings

% Data
dat_file = 'wber.csv';  % CSV file with prices by store, week, and UPC
                        % (download from https://www.chicagobooth.edu/research/kilts/datasets/dominicks under "Category Files")
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
V_hat_corr = sigma_hat.\V_hat./sigma_hat'; % Correlation matrix

% Results table
moment_tab = table;
moment_tab.estim = mu_hat';
moment_tab.se = sigma_hat;
moment_tab.corr = V_hat_corr(:,2:end);
moment_tab.Properties.RowNames = {'Freq', '2nd', '4th', 'Avg-abs'};
disp(moment_tab);

disp('Moments of price changes: std, kurt');
disp([sqrt(mu_hat(2)) sample_kurtosis]);



%% Estimation

disp(' ');
disp('ESTIMATION AND TESTING');
disp('Parameters: #products, vol, menu cost');

estim = estimate_model(mu_hat, V_hat);

% Results table
estim_tab = table;
estim_tab.justid = [estim.fullinfo.justid.theta';
                    estim.fullinfo.justid.se';
                    estim.indep.justid.theta';
                    estim.indep.justid.se';
                    estim.liminfo.justid.theta';
                    estim.liminfo.justid.se'];
estim_tab.overid_test = [estim.fullinfo.overid_test.errors(4)';
                    estim.fullinfo.overid_test.se(4)';
                    estim.indep.overid_test.errors(4)';
                    estim.indep.overid_test.se(4)';
                    estim.liminfo.overid_test.errors(4)';
                    estim.liminfo.overid_test.se(4)'];
estim_tab.eff =    [estim.fullinfo.eff.theta';
                    estim.fullinfo.eff.se';
                    estim.indep.eff.theta';
                    estim.indep.eff.se';
                    estim.liminfo.eff.theta';
                    estim.liminfo.eff.se'];
estim_tab.Properties.RowNames = {'Full-info estim', 'Full-info SE', 'Indep estim', 'Indep SE', 'Lim-info estim', 'Lim-info SE'};
disp(estim_tab);
                
% Display further results

disp(' ');
disp('Just-identified specification:');
disp('Ratio of SE: lim-info/full-info');
disp(estim.liminfo.justid.se'./estim.fullinfo.justid.se');
disp('Ratio of SE: lim-info/independence');
disp(estim.liminfo.justid.se'./estim.indep.justid.se');

disp(' ');
disp('Over-identification test:');
disp('Ratio of SE: lim-info/full-info');
disp(estim.liminfo.overid_test.se(4)'./estim.fullinfo.overid_test.se(4)');
disp('Ratio of SE: lim-info/independence');
disp(estim.liminfo.overid_test.se(4)'./estim.indep.overid_test.se(4)');
disp('Full-info p-value');
disp(2*normcdf(-abs(estim.fullinfo.overid_test.errors(4)./estim.fullinfo.overid_test.se(4))'));
disp('p-value under independence');
disp(2*normcdf(-abs(estim.indep.overid_test.errors(4)./estim.indep.overid_test.se(4))'));
disp('Worst-case p-value');
disp(2*normcdf(-abs(estim.liminfo.overid_test.errors(4)./estim.liminfo.overid_test.se(4))'));

disp(' ');
disp('Efficient specification:');
disp('Ratio of lim-info SE: just-ID/efficient');
disp(estim.liminfo.justid.se'./estim.liminfo.eff.se');
disp('Ratio of SE: lim-info-efficient/full-info-efficient');
disp(estim.liminfo.eff.se'./estim.fullinfo.eff.se');
disp('Ratio of SE: lim-info-efficient/independence-efficient');
disp(estim.liminfo.eff.se'./estim.indep.eff.se');


%% Simulation study, using empirically estimated parameters

if ~run_sim
    return;
end

disp(' ');
disp('SIMULATION STUDY');

% Ordering of procedures:
% 1: just-ID, full-info SE
% 2: just-ID, SE under independence
% 3: just-ID, lim-info SE
% 4: all moments, efficient, full-info SE
% 5: all moments, efficient, SE under independence
% 6: all moments, lim-info efficient, lim-info SE

% Preliminaries
theta_sim = estim.fullinfo.justid.theta; % True parameters in simulations
n_sim = n; % Sample size for simulations
[mu_sim, y_bar_sim] = moment_function(theta_sim); % True moments
quantiles_sim = quantile_price(theta_sim(1), y_bar_sim, quant_grid_sim); % Compute quantiles of price change distribution on grid (for quick simulation below)

% Run simulations
sim_theta_hat = nan(numrep_sim,3,6);
sim_se = sim_theta_hat;
sim_overid_tstat = nan(numrep_sim,3); % Last dimension: 1=efficient, 2=independent, 3=lim-info
sim_joint_pval = nan(numrep_sim,3);

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
        estim_sim = estimate_model(mu_hat_sim, V_hat_sim, theta_sim);
    catch ME
        fprintf('%s%d\n', 'Error for i=', i);
        warning(ME.msgtext);
        continue;
    end
    
    % Store estimates and SE
    sim_theta_hat(i,:,:) = [repmat(estim_sim.fullinfo.justid.theta,1,3) estim_sim.fullinfo.eff.theta estim_sim.indep.eff.theta estim_sim.liminfo.eff.theta];
    sim_se(i,:,:) = [estim_sim.fullinfo.justid.se estim_sim.indep.justid.se estim_sim.liminfo.justid.se estim_sim.fullinfo.eff.se estim_sim.indep.eff.se estim_sim.liminfo.eff.se];
    
    % Compute over-ID t-statistics
    sim_overid_tstat(i,:) = abs(estim_sim.fullinfo.overid_test.errors(4))./[estim_sim.fullinfo.overid_test.se(4) estim_sim.indep.overid_test.se(4) estim_sim.liminfo.overid_test.se(4)];
    
    % Joint test p-values
    sim_joint_pval(i,:) = [estim_sim.fullinfo.justid.joint_pval estim_sim.indep.justid.joint_pval estim_sim.liminfo.justid.joint_pval];
    
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
sim.joint_pval = sim_joint_pval;
sim.seeds = sim_seeds;
sim.elapsed_time = sim_elapsed_time;
clearvars sim_theta_hat sim_se sim_overid_tstat sim_joint_pval sim_seeds sim_elapsed_time;

% Coverage/rejection rates
cv_sim = norminv(1-alpha_sim/2); % Normal critical value
sim.theta_hat_error = sim.theta_hat-permute(theta_sim, [3 1 2]); % Estimation error
sim.cover = (abs(sim.theta_hat_error./sim.se)<cv_sim); % CI coverage indicator
sim.length = 2*sim.se*cv_sim; % CI length
sim.overid_reject = sim.overid_tstat>cv_sim; % Over-ID rejection indicator
sim.joint_reject = sim.joint_pval<alpha_sim; % Joint parameter test rejection indicator

sim.bias = permute(mean(sim.theta_hat_error,1,'omitnan'), [2 3 1]);
sim.var = permute(var(sim.theta_hat,0,1,'omitnan'), [2 3 1]);
sim.mse = sim.bias.^2+sim.var;
sim.cover_rate = permute(mean(sim.cover,1,'omitnan'), [2 3 1]);
sim.avg_length = permute(mean(sim.length,1,'omitnan'), [2 3 1]);
sim.med_length = permute(median(sim.length,1,'omitnan'), [2 3 1]);
sim.overid_reject_rate = mean(sim.overid_reject,1,'omitnan');
sim.joint_reject_rate = mean(sim.joint_reject,1,'omitnan');

% Results tables
row_names = {'Full-info', 'Indep', 'Lim-info'};

disp('RMSE relative to true parameters');
sim_tab.rmse = table;
sim_tab.rmse.justid = sqrt(sim.mse(:,1:3)')./theta_sim';
sim_tab.rmse.eff = sqrt(sim.mse(:,4:6)')./theta_sim';
sim_tab.rmse.Properties.RowNames = row_names;
disp(sim_tab.rmse);

disp('CI coverage rate');
sim_tab.cover = table;
sim_tab.cover.justid = sim.cover_rate(:,1:3)';
sim_tab.cover.eff = sim.cover_rate(:,4:6)';
sim_tab.cover.Properties.RowNames = row_names;
disp(sim_tab.cover);

disp('CI average length');
sim_tab.length = table;
sim_tab.length.justid = sim.avg_length(:,1:3)';
sim_tab.length.eff = sim.avg_length(:,4:6)';
sim_tab.length.Properties.RowNames = row_names;
disp(sim_tab.length);

disp('Over-ID rejection rate: avg-abs');
disp(sim.overid_reject_rate');

disp('Joint parameter test rejection rate:');
disp(sim.joint_reject_rate');

clearvars row_names;
