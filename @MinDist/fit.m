function res = fit(obj, varargin)

    % Minimum distance estimates and standard errors,
    % either with full-information moment var-cov matrix
    % or with limited-information individual moment variances
    
    
    %% Parse optional inputs

    ip = inputParser;
    addParameter(ip, 'transf', @(x) x, @(x) isa(x, 'function_handle'));
    addParameter(ip, 'weight_mat', [], @isnumeric);
    addParameter(ip, 'opt_init', [], @isnumeric);
    addParameter(ip, 'estim_fct', [], @(x) isa(x, 'function_handle') | isempty(x));
    addParameter(ip, 'eff', true,  @(x) islogical(x) | isnumeric(x));
    addParameter(ip, 'one_step', true,  @(x) islogical(x) | isnumeric(x));
    addParameter(ip, 'transf_deriv', [], @(x) isa(x, 'function_handle') | isempty(x));
    addParameter(ip, 'param_estim', [], @isnumeric);
    addParameter(ip, 'estim', [], @isnumeric);
    addParameter(ip, 'transf_jacob', [], @isnumeric);
    addParameter(ip, 'moment_fit', [], @isnumeric);
    addParameter(ip, 'moment_jacob', [], @isnumeric);
    parse(ip, varargin{:});
    
    
    %% Preliminaries
    
    % Check inputs
    assert(~isempty(ip.Results.param_estim) | ~isempty(ip.Results.opt_init) | ~isempty(ip.Results.estim_fct), ...
           'One of the following must be supplied: "param_estim", "opt_init", or "estim_fct"');
    
    % Transformation Jacobian function
    transf = ip.Results.transf;
    transf_deriv = deriv(ip.Results.transf_deriv, transf);
    
    % Determine weight matrix, if not supplied
    weight_mat = ip.Results.weight_mat;
    if obj.full_info && (ip.Results.eff || isempty(weight_mat))
        weight_mat = inv(obj.moment_varcov); % Full-info efficient weight matrix
    end
    if isempty(weight_mat)
        weight_mat = diag(1./diag(obj.moment_varcov)); % Ad hoc diagonal weight matrix
    end
    
    % Default estimation routine
    estim_fct = ip.Results.estim_fct;
    if isempty(estim_fct)
        opts = optimoptions('fminunc', 'Display', 'notify');
        estim_fct = @(W) fminunc(@(x) (obj.moment_estim-obj.moment_fct(x))'*W*(obj.moment_estim-obj.moment_fct(x)), ...
                                 ip.Results.opt_init, opts);
    end
    
    
    %%  Initial estimate of parameters, if not supplied
    
    param_estim = ip.Results.param_estim;
    if isempty(param_estim)
        param_estim = estim_fct(weight_mat);
    end
    
    % Transformation, moment function, and Jacobians at initial estimate
    [estim, transf_jacob, moment_fit, moment_jacob] ...
        = obj.estim_update(param_estim, transf, transf_deriv, ...
                           ip.Results.estim, ip.Results.transf_jacob, ip.Results.moment_fit, ip.Results.moment_jacob);
    moment_loadings = get_moment_loadings(moment_jacob, weight_mat, transf_jacob);
    estim_num = length(estim);
    
    
    %% Efficient estimates
    
    if ip.Results.eff % Update initial estimate to efficient estimate
        
        if obj.full_info % Full information
            
            if ip.Results.one_step % One-step estimation
                param_estim = obj.get_onestep(moment_fit, weight_mat, moment_jacob, param_estim);
                [estim, transf_jacob, moment_fit, moment_jacob] ...
                    = obj.estim_update(param_estim, transf, transf_deriv, ...
                                       [], [], [], []);
            else % Full optimization
                % Do nothing, since param_estim already contains estimates of interest
            end
            
        else % Limited information
            
            if estim_num>1 % If more than one parameter of interest, handle each separately by recursive call
                
                estim_init = estim;
                weight_mat_init = weight_mat;
                the_eye = eye(estim_num);
                
                estim = nan(estim_num,1);
                estim_se = nan(estim_num,1);
                moment_loadings = nan(obj.moment_num,estim_num);
                weight_mat = cell(1,estim_num);
                
                for i=1:estim_num % Loop over each parameter of interest
                    the_res = obj.fit('transf', @(x) the_eye(m,:)*transf(x), ...
                                      'weight_mat', weight_mat_init, ...
                                      'estim_fct', estim_fct, ...
                                      'eff', true, ...
                                      'one_step', ip.Results.one_step, ...
                                      'param_estim', param_estim, ...
                                      'estim', estim_init(i), ...
                                      'transf_jacob', transf_jacob(i,:), ...
                                      'moment_fit', moment_fit, ...
                                      'moment_jacob', moment_jacob);
                    estim(i) = the_res.estim;
                    estim_se(i) = the_res.estim_se;
                    moment_loadings(:,i) = the_res.moment_loadings;
                    weight_mat{i} = the_res.weight_mat;
                end
                
            else % If only single parameter of interest
                
                [estim_se, moment_loadings, weight_mat] ...
                    = obj.worstcase_eff(moment_jacob, transf_jacob, weight_mat);
                if ip.Results.one_step % One-step estimation
                    estim = obj.get_onestep(moment_fit, [], moment_loadings, estim);
                else % Full optimization estimation
                    param_estim = estim_weight(weight_mat);
                    estim = transf(param_estim);
                end
            end
            
        end
        
    end
    
    
    %% Start building results structure
    
    res = struct;
    res.estim = estim;
	res.param_estim = param_estim;
	res.weight_mat = weight_mat;
    res.moment_fit = moment_fit;
	res.moment_jacob = moment_jacob;
    res.moment_loadings = moment_loadings;
	res.transf_jacob = transf_jacob;
    res.eff = ip.Results.eff;
	res.estim_num = estim_num;
    
    
    %% Standard errors
    
    if obj.full_info % Full information
        
        estim_varcov = moment_loadings' * obj.moment_varcov * moment_loadings;
        estim_se = sqrt(diag(estim_varcov));
        res.estim_varcov = estim_varcov;
        
    else % Limited information
        
        if ip.Results.eff % SE have already been computed above for worst-case efficient estimates
            % Do nothing, since standard errors have already been computed above
        else
            [estim_se, worstcase_varcov] = obj.worstcase_se(moment_loadings);
            res.worstcase_varcov = worstcase_varcov;
        end
        
    end
    
    res.estim_se = estim_se;
    
end