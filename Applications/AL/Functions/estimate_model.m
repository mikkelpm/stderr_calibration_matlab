function estim = estimate_model(mu_hat, V_hat)

    % Estimate and test Alvarez & Lippi (ECMA 2014) model in various ways

    
    estim = struct;
    sigma_hat = sqrt(diag(V_hat)); % Marginal standard errors of moments
    p = length(mu_hat);
    

    %% Estimate model: just-identified moments (first three)
    
    % Full information
    W_justid = blkdiag(eye(3),zeros(p-3,p-3)); % Weight matrix for just-identified estimation
    res_justid_fullinfo = WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                      'W', W_justid, 'V', V_hat, 'eff', false);
    estim.justid.theta_hat = res_justid_fullinfo.theta; % Estimates
    estim.justid.fullinfo_se = res_justid_fullinfo.r_theta_se; % SE

    % Assuming independence
    res_justid_indep = WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                   'W', W_justid, 'V', diag(diag(V_hat)), 'eff', false);
    estim.justid.indep_se = res_justid_indep.r_theta_se; % SE
    
    % Worst case
    res_justid_wc = WorstCaseSE(@moment_function, mu_hat, sigma_hat, @(W) param_closed_form(mu_hat(1:3)), ...
                                'W', W_justid, 'eff', false);
    estim.justid.wc_se = res_justid_wc.r_theta_se; % SE


    %% Test of over-identifying restrictions (fourth moment)

    % Full information
    estim.overid_test.errors = res_justid_fullinfo.overid.errors; % Errors in matching moments
    estim.overid_test.fullinfo_se = res_justid_fullinfo.overid.se; % SE for moment errors
    
    % Assuming independence
    estim.overid_test.indep_se = res_justid_indep.overid.se; % SE for moment errors
    
    % Worst case
    estim.overid_test.wc_se = res_justid_wc.overid.se; % SE for moment errors
    
    % Can't test validity of first three moments (used for just-ID estimation)
    estim.overid_test.fullinfo_se(1:3) = nan;
    estim.overid_test.indep_se(1:3) = nan;
    estim.overid_test.wc_se(1:3) = nan;

    
    %% Estimation by worst-case-efficient moment selection

    res_eff_wc = WorstCaseSE(@moment_function, mu_hat, sigma_hat, @(W) param_closed_form(mu_hat(1:3)), ...
                             'W', W_justid, 'eff', true, 'one_step', true);
    estim.wceff.theta_hat = res_eff_wc.r_theta; % Estimate
    estim.wceff.se = res_eff_wc.r_theta_se; % SE
    estim.wceff.x_hat = res_eff_wc.x_hat; % Moment loadings


    %% Full-information efficient estimation

    res_eff_fullinfo = WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                   'V', V_hat, 'eff', true, 'one_step', true);
    estim.fullinfo.theta_hat = res_eff_fullinfo.theta; % Estimate
    estim.fullinfo.se = res_eff_fullinfo.r_theta_se; % SE
    
    
end