clear all;
addpath('Functions');
addpath('../../Main');
addpath('../../Supporting');

% Application based on Alvarez & Lippi (ECMA 2014)


%% Settings

dat_file = 'wber.csv';  % CSV file with prices by store, week, and UPC
select_store = 122;     % Single store ID to keep for the data analysis
drop_quantile = 0.99;   % Drop absolute price changes above this quantile


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
sample_kurtosis = mu_hat(3)/mu_hat(2)^2;

disp('Fraction of nonzero price changes per week');
disp(mu_hat(1));

disp('Moments of price changes: std, kurt');
disp([sqrt(mu_hat(2)) sample_kurtosis]);

disp('Moments of absolute price changes: mean, std');
disp([mu_hat(4) sqrt(mu_hat(5))]);


%% Estimation

disp(' ');
disp('ESTIMATION');

estim = estimate_parameters(mu_hat, V_hat);

disp(' ');
disp('Just-identified estimates: #prod, vol, sqrt(menu cost)');
disp(estim.justid.theta_hat');
disp('Actual SE');
disp(estim.justid.se');
disp('SE under independence');
disp(estim.justid.indep_se');
disp('Worst-case SE');
disp(estim.justid.wcse');
disp('Ratio of SE: actual/worst-case');
disp(estim.justid.se'./estim.justid.wcse');
disp('Ratio of SE: independence/worst-case');
disp(estim.justid.indep_se'./estim.justid.wcse');

disp(' ');
disp('Over-ID test moment errors: avg#changes, std, 4th, avg-abs, std-abs');
disp(estim.overid_test.moment_errors');
disp('Actual SE');
disp(estim.overid_test.se');
disp('Worst-case SE');
disp(estim.overid_test.wcse');
disp('Ratio of SE: actual/worst-case');
disp(estim.overid_test.se'./estim.overid_test.wcse');
disp('Actual p-values');
disp(2*(1-normcdf(abs(estim.overid_test.moment_errors./estim.overid_test.se)')));
disp('Worst-case p-values');
disp(2*(1-normcdf(abs(estim.overid_test.moment_errors./estim.overid_test.wcse)')));

disp(' ');
disp('Worst-case optimal estimates: #prod, vol, sqrt(menu cost)');
disp(estim.wcopt.theta_hat');
disp('Worst-case SE');
disp(estim.wcopt.wcse');
disp('Ratio of worst-case SE: optimal/just-ID');
disp(estim.wcopt.wcse'./estim.justid.wcse');

disp(' ');
disp('Full-information efficient estimates: #prod, vol, sqrt(menu cost)');
disp(estim.fullinfo.theta_hat');
disp('SE');
disp(estim.fullinfo.se');
disp('Ratio of SE: full-info/worst-case-optimal');
disp(estim.fullinfo.se'./estim.wcopt.wcse');

