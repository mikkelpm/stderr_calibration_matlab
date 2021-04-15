function estim = estimate_model(mu_hat, V_hat)

    % Estimate and test Alvarez & Lippi (ECMA 2014) model in various ways

    
    estim = struct;
    sigma_hat = sqrt(diag(V_hat)); % Marginal standard errors of moments
    

    %% Estimate model: just-identified moments

    % Just-identified estimates, using only three moments
    estim.justid.theta_hat = param_closed_form(mu_hat(1:3)); % Closed-form expression

    % Standard errors, actual, under independence, and worst case
    jacob_justid = FiniteDiff(@(x) param_closed_form(x), mu_hat(1:3)'); % Jacobian in closed-form expression for hat{theta}
    estim.justid.se = sqrt(diag(jacob_justid*V_hat(1:3,1:3)*jacob_justid')); % Actual SE
    estim.justid.indep_se = sqrt(diag(jacob_justid*diag(sigma_hat(1:3).^2)*jacob_justid')); % SE if independent
    estim.justid.wcse = abs(jacob_justid)*sigma_hat(1:3); % Worst-case SE


    %% Test over-identifying restrictions

    p = length(mu_hat);
    h = @(theta) moment_function(theta);
    G_fct = @(theta) FiniteDiff(h, theta);

    h_theta_hat_justid = h(estim.justid.theta_hat); % Moment function at just-ID estimates
    estim.overid_test.moment_errors = mu_hat'-h_theta_hat_justid; % Moment errors
    G_hat_justid = G_fct(estim.justid.theta_hat); % hat{G}
    M_hat_justid = eye(p) - G_hat_justid*[jacob_justid, zeros(3,p-3)]; % Jacobian for transforming sample moments into moment errors
    estim.overid_test.se = sqrt(diag(M_hat_justid*V_hat*M_hat_justid')); % SE for moment errors, actual
    estim.overid_test.indep_se = sqrt(diag(M_hat_justid*diag(sigma_hat.^2)*M_hat_justid')); % SE if independent
    estim.overid_test.wcse = abs(M_hat_justid)*sigma_hat; % Worst-case SE
    
    % Can't test validity of first three moments (used for just-ID estimation)
    estim.overid_test.se(1:3) = nan;
    estim.overid_test.indep_se(1:3) = nan;
    estim.overid_test.wcse(1:3) = nan;

    
    %% Estimation by worst-case-optimal moment selection

    k = length(estim.justid.theta_hat);

    estim.wcopt.theta_hat = nan(k,1);
    estim.wcopt.wcse = nan(k,1);
    estim.wcopt.x_hat = nan(k,p);
    the_eye = eye(k);

    for i=1:k % Loop over parameters
        [~, res_onestep] ...
            = HTest_Single_WorstCaseOptW_Step2(V_hat, ...
                G_hat_justid, G_fct, the_eye(:,i), 1, ...
                false, h, [], ...
                true, mu_hat', estim.justid.theta_hat, h_theta_hat_justid); % Worst-case optimal estimates (one-step update from just-ID estimates)
        estim.wcopt.theta_hat(i) = res_onestep.lcomb; % Parameters estimates
        estim.wcopt.wcse(i) = res_onestep.stderr_check; % SE
        estim.wcopt.x_hat(i,:) = res_onestep.x_check; % Loadings on moments
    end


    %% Full-information efficient estimation

    estim.fullinfo.theta_hat = getOnestep(h_theta_hat_justid, inv(V_hat), G_hat_justid, estim.justid.theta_hat, mu_hat'); % Efficient estimates (one-step update from just-ID estimates)
    estim.fullinfo.se = sqrt(diag(inv(G_hat_justid'*(V_hat\G_hat_justid)))); % SE
    
    
end