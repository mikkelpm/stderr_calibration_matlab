function estim = estimate_model(mu_hat, V_hat, varargin)

    % Estimate and test Alvarez & Lippi (ECMA 2014) model in various ways
    
    
    %% Specifications
    
    obj = cell(1,3);
    obj_names = {'fullinfo', 'indep', 'liminfo'};
    
    obj{1} = MinDist(@moment_function, mu_hat, 'moment_varcov', V_hat); % Full information
    obj{2} = MinDist(@moment_function, mu_hat, 'moment_varcov', diag(diag(V_hat))); % Assuming independence
    obj{3} = MinDist(@moment_function, mu_hat, 'moment_se', sqrt(diag(V_hat))); % Limited information
    
    
    %% Estimate and test each specification
    
    estim = struct;
    p = length(mu_hat);
    W_justid = blkdiag(eye(3),zeros(p-3,p-3)); % Weight matrix for just-identified estimation
    
    for j=1:length(obj)

        the_estim = struct;
        
        % Just-identified specification (use first three moments only)
        res_justid = obj{j}.fit('estim_fct', @(W) param_closed_form(mu_hat(1:3)), ...
                                'weight_mat', W_justid, 'eff', false);
        the_estim.justid.theta = res_justid.estim; % Estimates
        the_estim.justid.se = res_justid.estim_se; % SE

        % Test of over-identifying restrictions (fourth moment)
        res_overid = obj{j}.overid(res_justid, 'joint', false);
        the_estim.overid_test.errors = res_overid.moment_error; % Errors in matching moments
        the_estim.overid_test.se = res_overid.moment_error_se; % SE for moment errors
        the_estim.overid_test.se(1:3) = nan; % Can't test validity of first three moments (used for just-ID estimation)

        % Estimate model: efficient use of moments
        res_eff = obj{j}.fit('estim_fct', @(W) param_closed_form(mu_hat(1:3)), ...
                             'weight_mat', W_justid, 'eff', true, 'one_step', true);
        the_estim.eff.theta = res_eff.estim; % Estimate
        the_estim.eff.se = res_eff.estim_se; % SE
        the_estim.eff.moment_loadings = res_eff.moment_loadings; % SE

        % Joint hypothesis test based on just-identified moments
        if ~isempty(varargin) % If parameter values are supplied...
            theta0 = varargin{1}(:); % Hypothesized parameter values
            res_justid_restr = obj{j}.fit('transf', @(x) theta0-x, 'estim_fct', @(W) param_closed_form(mu_hat(1:3)), ...
                                          'weight_mat', W_justid, 'eff', false);
            res_justid_test = obj{j}.test(res_justid_restr, 'joint', true);
            the_estim.justid.joint_pval = res_justid_test.joint_pval;
        end
        
        estim.(obj_names{j}) = the_estim;
        
    end
    
    
end