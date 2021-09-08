function res = test(obj, estim_res, varargin)

    % Test whether transformed parameters equal zero

    
    %% Parse optional inputs

    ip = inputParser;
    addRequired(ip, 'estim_res', @isstruct);
    addParameter(ip, 'joint', true, @islogical);
    addParameter(ip, 'test_weight_mat', [], @isnumeric);
    parse(ip, estim_res, varargin{:});

    
    %% t-statistics
    
    tstat = estim_res.estim./estim_res.estim_se;
    tstat_pval = 2*normcdf(-abs(tstat));
    
    res = struct;
    res.tstat = tstat;
    res.tstat_pval = tstat_pval;
    
    if ~ip.Results.joint
        return;
    end
    
    
    %% Joint test
    
    % Weight matrix for joint test statistics
    test_weight_mat = ip.Results.test_weight_mat;
    if isempty(test_weight_mat)
        if obj.full_info % Full information
            test_weight_mat = inv(estim_res.estim_varcov);
        else % Limited information
            test_weight_mat = inv(estim_res.moment_loadings'*diag(diag(obj.moment_varcov))*estim_res.moment_loadings);
        end
    end
    
    % Check dimensions
    assert(size(test_weight_mat,1)==estim_res.estim_num & size(test_weight_mat,2)==estim_res.estim_num, ...
           'Dimension of "test_weight_mat" is wrong');
    
    % Test statistic
    joint_stat = estim_res.estim' * test_weight_mat * estim_res.estim;
    
    % p-value
    if obj.full_info % Full information
        joint_pval = 1-chi2cdf(joint_stat, estim_res.estim_num);
    else % Limited information
        [max_trace, max_trace_varcov] = obj.solve_sdp(estim_res.moment_loadings*test_weight_mat*estim_res.moment_loadings');
        joint_pval = 1-chi2cdf(joint_stat/max_trace, 1);
        if joint_pval>0.215 % Test can only be used at significance levels < 0.215
            joint_pval = 1;
        end
        res.max_trace = max_trace;
        res.max_trace_varcov = max_trace_varcov;
    end
    
    res.test_weight_mat = test_weight_mat;
    res.joint_stat = joint_stat;
    res.joint_pval = joint_pval;
    
        

end