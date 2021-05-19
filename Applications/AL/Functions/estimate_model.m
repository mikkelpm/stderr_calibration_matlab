function estim = estimate_model(mu_hat, V_hat, varargin)

    % Estimate and test Alvarez & Lippi (ECMA 2014) model in various ways

    
    estim = struct;
    sigma_hat = sqrt(diag(V_hat)); % Marginal standard errors of moments
    p = length(mu_hat);
    

    %% Estimate model: just-identified moments (first three)
    
    % Full information
    W_justid = blkdiag(eye(3),zeros(p-3,p-3)); % Weight matrix for just-identified estimation
    res_justid_fullinfo = WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                      'W', W_justid, 'V', V_hat, 'eff', false);
    estim.justid.theta = res_justid_fullinfo.theta; % Estimates
    estim.justid.fullinfo.se = res_justid_fullinfo.r_theta_se; % SE

    % Assuming independence
    res_justid_indep = WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                   'W', W_justid, 'V', diag(diag(V_hat)), 'eff', false);
    estim.justid.indep.se = res_justid_indep.r_theta_se; % SE
    
    % Worst case
    res_justid_wc = WorstCaseSE(@moment_function, mu_hat, sigma_hat, @(W) param_closed_form(mu_hat(1:3)), ...
                                'W', W_justid, 'eff', false);
    estim.justid.wc.se = res_justid_wc.r_theta_se; % SE


    %% Test of over-identifying restrictions (fourth moment)

    % Full information
    estim.overid_test.errors = res_justid_fullinfo.overid.errors; % Errors in matching moments
    estim.overid_test.fullinfo.se = res_justid_fullinfo.overid.se; % SE for moment errors
    
    % Assuming independence
    estim.overid_test.indep.se = res_justid_indep.overid.se; % SE for moment errors
    
    % Worst case
    estim.overid_test.wc.se = res_justid_wc.overid.se; % SE for moment errors
    
    % Can't test validity of first three moments (used for just-ID estimation)
    estim.overid_test.fullinfo.se(1:3) = nan;
    estim.overid_test.indep.se(1:3) = nan;
    estim.overid_test.wc.se(1:3) = nan;

    
    %% Estimate model: efficient use of moments

    % Full information
    res_eff_fullinfo = WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                   'V', V_hat, 'eff', true, 'one_step', true);
    estim.eff.fullinfo.theta = res_eff_fullinfo.theta; % Estimate
    estim.eff.fullinfo.se = res_eff_fullinfo.r_theta_se; % SE
    
    % Assuming independence
    res_eff_indep = WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                   'V', diag(diag(V_hat)), 'eff', true, 'one_step', true);
    estim.eff.indep.theta = res_eff_indep.theta; % Estimate
    estim.eff.indep.se = res_eff_indep.r_theta_se; % SE
    
    % Worst case
    res_eff_wc = WorstCaseSE(@moment_function, mu_hat, sigma_hat, @(W) param_closed_form(mu_hat(1:3)), ...
                             'W', W_justid, 'eff', true, 'one_step', true);
    estim.eff.wc.theta = res_eff_wc.r_theta; % Estimate
    estim.eff.wc.se = res_eff_wc.r_theta_se; % SE
    estim.eff.wc.x_hat = res_eff_wc.x_hat; % Moment loadings
    
    
    %% Joint hypothesis test based on just-identified moments
    
    if ~isempty(varargin) % If parameter values are supplied...
        
        theta0 = varargin{1}(:); % Hypothesized parameter values
        
        % Full information
        res_justid_fullinfo_joint = WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                                'r', @(x)x-theta0, 'W', W_justid, 'V', V_hat, 'eff', false, 'joint', true);
        estim.justid.fullinfo.joint_pval = res_justid_fullinfo_joint.joint.p_val;
        
        % Assuming independence
        res_justid_indep_joint =    WorstCaseSE(@moment_function, mu_hat, [], @(W) param_closed_form(mu_hat(1:3)), ...
                                                'r', @(x)x-theta0, 'W', W_justid, 'V', diag(diag(V_hat)), 'eff', false, 'joint', true);
        estim.justid.indep.joint_pval = res_justid_indep_joint.joint.p_val;
        
        % Worst case
        res_justid_wc_joint =       WorstCaseSE(@moment_function, mu_hat, sigma_hat, @(W) param_closed_form(mu_hat(1:3)), ...
                                                'r', @(x)x-theta0, 'W', W_justid, 'eff', false, 'joint', true);
        estim.justid.wc.joint_pval = res_justid_wc_joint.joint.p_val;
        
    end
    
    
end