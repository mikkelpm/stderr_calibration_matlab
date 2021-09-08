function res = overid(obj, estim_res, varargin)

    % Over-identification test
    
    
    %% Parse optional inputs

    ip = inputParser;
    addRequired(ip, 'estim_res', @isstruct);
    addParameter(ip, 'joint', true, @islogical);
    parse(ip, estim_res, varargin{:});
    
    assert(isnumeric(estim_res.weight_mat), 'Estimation results must be based on a single weight matrix');
    
    
    %% Moment errors and their standard errors
    
    % Errors in fitting moments
    moment_error = obj.moment_estim - estim_res.moment_fit;
    
    % Standard errors for moment errors
    M = eye(obj.moment_num) ...
        - get_moment_loadings(estim_res.moment_jacob, estim_res.weight_mat, estim_res.moment_jacob)';
    
    the_estim_res = obj.fit('transf', @(x) x, ...
                            'eff', false, ...
                            'weight_mat', eye(obj.moment_num), ...
                            'param_estim', estim_res.param_estim, ...
                            'estim', moment_error, ...
                            'transf_jacob', M, ...
                            'moment_fit', obj.moment_estim, ...
                            'moment_jacob', eye(obj.moment_num));
    % Only the inputs "weight_mat", "transf_jacob", and "moment_jacob" are
    % actually used to calculate the standard errors - the other inputs are
    % only provided to avoid unnecessary computations
    
    
    %% Test statistic and p-value
    
    the_estim_res.eff = true; % Trick test() function below into not throwing an error in full-info case
    the_test_res = obj.test(the_estim_res, 'joint', ip.Results.joint, 'test_weight_mat', estim_res.weight_mat);
    
    res = struct;
    res.moment_error = moment_error;
    res.moment_error_se = the_estim_res.estim_se;
    res.tstat = the_test_res.tstat;
    res.tstat_pval = the_test_res.tstat_pval;
    
    if ip.Results.joint
        res.joint_stat = the_test_res.joint_stat;
        res.joint_pval = the_test_res.joint_pval;
        if obj.full_info % Adjust degrees of freedom
            assert(estim_res.eff, 'Full-information joint test requires using the efficient weight matrix');
            res.joint_pval = 1-chi2cdf(the_test_res.joint_stat, obj.moment_num-length(estim_res.param_estim));
        else
            res.max_trace = the_test_res.max_trace;
            res.max_trace_varcov = the_test_res.max_trace_varcov;
        end
    end
    

end