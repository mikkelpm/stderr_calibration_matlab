function res = WorstCaseSE(h, mu, theta_estim_fct, varargin)

    % Worst-case standard errors for minimum distance estimates
    % without knowledge of the correlation matrix for the matched moments
    
    % If desired, also computes:
    % - worst-case optimal estimates
    % - full-information efficient estimates
    % - over-identification test for each individual moments
    
    % Reference:
    % Cocci, Matt & Mikkel Plagborg-Moller, "Standard Errors for Calibrated
    % Parameters", https://scholar.princeton.edu/mikkelpm/calibration

    
    %% Parse inputs

    p = length(mu);
    
    ip = inputParser;
    
    % Required inputs
    addRequired(ip, 'h', @(x) isa(x, 'function_handle'));
        % Function mapping parameters theta (k x 1) into moments mu (p x 1)
    addRequired(ip, 'mu', @isnumeric);
        % Estimated moments mu (p x 1)
    addRequired(ip, 'theta_estim_fct', @(x) isa(x, 'function_handle') | isempty(x));
        % Function mapping mu (p x 1) and weight matrix W (p x p) into minimum distance estimate theta (k x 1)
    
    % Optional inputs
    addParameter(ip, 'lambda', [], @isnumeric);
        % Linear combinations (k x r) of theta of interest, each column represents a different linear combination of length k
        % Default: lambda=eye(k), i.e., interested in each element of theta separately
    addParameter(ip, 'sigma', ones(p,1), @isnumeric);
        % Standard errors of mu estimate (p x 1), only needed for worst-case analysis
        % Default: vector of ones
    addParameter(ip, 'W', [], @isnumeric);
        % Minimum distance weight matrix (p x p)
        % Default: diagonal matrix with entries sigma_j^{-2}
    addParameter(ip, 'V', [], @isnumeric);
        % Full-information var-cov matrix (p x p) of mu - if supplied, yields full-information analysis
        % Default: empty, yielding worst-case analysis
    addParameter(ip, 'theta', [], @isnumeric);
        % Initial consistent estimate of theta (k x 1) - if supplied, over-rides initial estimation step
        % Default: empty
    addParameter(ip, 'opt', true,  @(x) islogical(x) | isnumeric(x));
        % true: optimal estimation (either full-information-efficient or worst-case-optimal)
        % Default: true
    addParameter(ip, 'one_step', true,  @(x) islogical(x) | isnumeric(x));
        % true: one-step optimal estimate, false: full re-optimization via "theta_estim_fct" function
        % Default: true
    addParameter(ip, 'overid', true, @(x) islogical(x) | isnumeric(x));
        % true: compute over-identification tests
        % Default: true
    addParameter(ip, 'zero_thresh', 1e-8, @(x) isnumeric(x) && isscalar(x));
        % Numerical threshold for determining whether eigenvalues are zero
        % Default: 1e-8
    
    parse(ip, h, mu, theta_estim_fct, varargin{:}); % Parse inputs
    
    
    %% Check inputs
    
    res = struct;
    
    res.mu = mu(:);
    res.sigma = ip.Results.sigma(:);
    res.W = ip.Results.W;
    res.V = ip.Results.V;
    
    if isempty(res.W)
        res.W = diag(1./res.sigma.^2); % Default weight matrix
    end
    
    assert(length(res.sigma)==p, 'Dimension of sigma is incorrect');
    assert(isequal(size(res.W),[p,p]), 'Dimensions of W are incorrect');
    assert(isempty(res.V) || isequal(size(res.V),[p,p]), 'Dimensions of V are incorrect');
    assert(ip.Results.zero_thresh>=0, 'zero_thresh must be non-negative');
    
    res.full_info = ~isempty(res.V); % Full-information estimate? (Requires matrix V.)
    res.opt = ip.Results.opt; % Optimal estimate?
    res.one_step = ip.Results.one_step; % One-step optimal estimate?
    
    
    %% Initial estimate
    
    if ~isempty(ip.Results.theta) % If initial estimate is already supplied
        res.theta = ip.Results.theta;
    else % Otherwise, compute initial estimate given mu and W
        res.theta = theta_estim_fct(res.mu, res.W);
    end
    res.theta = res.theta(:);
    res.h_theta = h(res.theta);
    k = length(res.theta);
    res.G = ComputeG(h, res.theta); % Moment function Jacobian at initial estimate
    assert(isequal(size(res.G),[p,k]), 'Jacobian of moment function has incorrect dimensions');
    
    % Check/define lambda
    res.lambda = ip.Results.lambda;
    if isempty(res.lambda)
        res.lambda = eye(k); % Default: Compute SE for each element of theta
    elseif size(res.lambda,1)==1
        res.lambda = res.lambda';
    end
    assert(size(res.lambda,1)==k, 'Dimensions of lambda are incorrect');
    
    res.lambda_theta = res.lambda'*res.theta; % Linear combinations of interest, initial estimate
    
    
    %% Updated estimate
    
    if res.opt % Update initial estimate to optimal estimate
        
        % Store initial estimates
        res.theta_init = res.theta;
        res.h_theta_init = res.h_theta;
        res = rmfield(res, {'theta', 'h_theta'});
        
        if res.full_info % Full information analysis
            
            res.W = inv(res.V); % Efficient weight matrix
            
            if res.one_step % One-step update
                res.theta = getOnestep(res.h_theta_init, res.W, res.G, res.theta_init, res.mu);
            else % Full optimization update
                res.theta = theta_estim_fct(res.mu, res.W);
            end
            
            res.lambda_theta = res.lambda'*res.theta; % Linear combinations of interest
            res.G = ComputeG(h, res.theta); % Update Jacobian with new estimate
            res.h_theta = h(res.theta); % Update h(hat{theta})
            
        else % Worst case analysis
            
            nl = size(res.lambda,2);
            
            if nl>1 % If there is more than one parameter of interest
                
                res.lambda_theta = nan(nl,1);
                res.lambda_theta_se = nan(nl,1);
                res.x_hat = nan(nl,p);
                
                for il=1:nl % Loop over each parameter of interest
                    % Recursive call to function with single vector lambda
                    the_res = WorstCaseSE(h, res.mu, theta_estim_fct, ...
                                          'lambda', res.lambda(:,il), ...
                                          'sigma', res.sigma, 'W', res.W, ...
                                          'theta', res.theta_init, ...
                                          'opt', true, 'one_step', res.one_step, ...
                                          'overid', false, 'zero_thresh', ip.Results.zero_thresh);
                    res.lambda_theta(il) = the_res.lambda_theta;
                    res.lambda_theta_se(il) = the_res.lambda_theta_se;
                    res.x_hat(il,:) = the_res.x_hat;
                end
                
            else % Otherwise, compute worst-case optimum for the single parameter of interest
                
                % Worst-case optimal weighting and SE
                [res.x_hat, res.lambda_theta_se] = ComputeWorstCaseOpt_Single(res.sigma, res.G, res.lambda, ip.Results.zero_thresh);
                
                % Update point estimate
                if res.one_step % One-step update
                    res.lambda_theta = getOnestep(res.h_theta_init, [], res.x_hat, res.lambda'*res.theta_init, res.mu);
                else % Full optimization update
                    % Weight matrix only puts weight on k moments with largest loadings in hat{x}
                    [~,I]=sort(abs(x_hat)); % Sort elements in hat{x} by magnitude
                    res.W(I(1:p-k),:) = 0;
                    res.W(:,I(1:p-k)) = 0;
                    res.theta = theta_estim_fct(res.mu, res.W); % Update estimate
                    res.h_theta = h(res.theta); % Update h(hat{theta})
                    res.lambda_theta = res.lambda'*res.theta; % Linear combination of interest
                end
                
            end
            
        end
        
    end
    
    
    %% SE
    
    if res.full_info % Full information analysis
        
        res.lambda_theta_varcov = ComputeCMDAvar(res.G,res.W,res.V,res.lambda'); % Var-cov matrix
        res.lambda_theta_se = sqrt(diag(res.lambda_theta_varcov)); % SE
        
    else % Worst case analysis
        
        if ~res.opt % SE have already been computed above for worst-case optimal estimates
            
            res.x_hat = res.W*res.G*((res.G'*res.W*res.G)\res.lambda); % Moment loadings
            res.lambda_theta_se = abs(res.x_hat)'*res.sigma; % SE
            
        end
        
    end
    
    
    %% Over-identification test
    
    if ip.Results.overid
        
        res.overid = struct;
        
        % Errors in matching moments
        if isfield(res, 'h_theta')
            res.overid.errors = res.mu-res.h_theta;
        else
            res.overid.errors = res.mu-res.h_theta_init;
        end
        
        % SE for moment errors
        M = eye(p) - res.W*res.G*((res.G'*res.W*res.G)\res.G');
        the_res = WorstCaseSE(@(x) x, res.mu, [], ...
                              'lambda', M, ...
                              'sigma', res.sigma, 'W', eye(p), ...
                              'V', res.V, 'theta', res.mu, ...
                              'opt', false, 'overid', false, ...
                              'zero_thresh', ip.Results.zero_thresh); % Compute standard errors for M'*hat{mu}
        res.overid.se = the_res.lambda_theta_se;
        
        % Test statistics and p-values
        res.overid.t_stats = res.overid.errors./res.overid.se;
        res.overid.p_vals = 2*normcdf(-abs(res.overid.t_stats));
        
    end
    
end