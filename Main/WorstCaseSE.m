function res = WorstCaseSE(h, mu, sigma, theta_estim_fct, varargin)

    % Worst-case standard errors for minimum distance estimates
    % without knowledge of the correlation matrix for the matched moments
    
    % If desired, also computes:
    % - worst-case optimal estimates
    % - full-information efficient estimates
    % - over-identification test for each individual moment
    
    % Reference:
    % Cocci, Matthew D. & Mikkel Plagborg-Moller, "Standard Errors for Calibrated Parameters"
    % https://scholar.princeton.edu/mikkelpm/calibration

    
    % Inputs: see below
    
    % Outputs: structure with the following fields
    % (some fields may be missing in special cases)
    % - theta:                  point estimates hat{theta} (k x 1)
    % - h_theta:                moment function evaluated at hat{theta} (p x 1)
    % - r_theta:                point estimates r(hat{theta}) (m x 1)
    % - r_theta_se:             SE of r(hat{theta}) (m x 1)
    % - r_theta_varcov:         var-cov matrix of r(hat{theta}) (m x m), if full info
    % - G:                      Jacobian of h(theta) (p x k)
    % - R:                      Jacobian of r(theta) (m x k)
    % - x_hat:                  loadings of r(hat{theta}) on moments hat{mu} (p x m)
    % - overid:                 structure for over-identification tests for individual moments
    %                           with the following fields
    % -- errors:                moment errors hat{mu}-h(hat{theta})
    % -- se:                    SE of moment errors
    % -- t_stats:               t-statistics for testing moment errors equal zero
    % -- p_vals:                p-values for t-tests
    
    
    %% Parse inputs

    p = length(mu);
    
    ip = inputParser;
    
    % Required inputs
    addRequired(ip, 'h', @(x) isa(x, 'function_handle'));
        % Function mapping parameters theta (k x 1) into moments mu (p x 1)
    addRequired(ip, 'mu', @isnumeric);
        % Estimated moments hat{mu} (p x 1)
    addRequired(ip, 'sigma', @isnumeric);
        % Standard errors of hat{mu} (p x 1), may be empty if "V" is supplied (see below)
    addRequired(ip, 'theta_estim_fct', @(x) isa(x, 'function_handle') | isempty(x));
        % Function mapping weight matrix W (p x p) into minimum distance estimate theta (k x 1),
        % e.g., theta_estim_fct = @(W) fminunc(@(theta) (mu-h(theta))'*W*(mu-h(theta)), x0)
    
    % Optional inputs
    addParameter(ip, 'r', @(x) x, @(x) isa(x, 'function_handle') | isempty(x));
        % Function that returns transformed parameters r(theta) of interest, given theta, may be vector-valued (m x 1)
        % Default: r(theta)=theta
    addParameter(ip, 'W', diag(1./sigma.^2), @isnumeric);
        % Minimum distance weight matrix (p x p)
        % Default: diagonal matrix with entries sigma_j^{-2}
    addParameter(ip, 'V', [], @isnumeric);
        % Full-information var-cov matrix (p x p) of mu - if supplied, yields full-information analysis
        % Default: empty, yielding worst-case analysis
    addParameter(ip, 'theta', [], @isnumeric);
        % Initial consistent estimate of theta (k x 1) - if supplied, over-rides initial estimation step
        % Default: empty
    addParameter(ip, 'G_fct', @(x) ComputeG(h,x), @(x) isa(x, 'function_handle') || isnumeric(x));
        % Function that returns Jacobian (p x k) of h(.) with respect to theta, given theta
        % Default: numerically differentiate h(.)
    addParameter(ip, 'R_fct', [], @(x) isa(x, 'function_handle') || isnumeric(x));
        % Function that returns Jacobian (m x k) of r(.) with respect to theta, given theta
        % Default: numerically differentiate r(.)
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
    
    parse(ip, h, mu, sigma, theta_estim_fct, varargin{:}); % Parse inputs
    
    
    %% Check inputs
    
    res = struct;
    
    res.mu = mu(:);
    res.sigma = sigma(:);
    res.W = ip.Results.W;
    res.V = ip.Results.V;
    
    if isempty(res.W) && ~isempty(res.V)
        res.W = inv(res.V); % Full-info efficient weight matrix
    end
    
    G_fct = ip.Results.G_fct;
    if isnumeric(G_fct)
        G_fct = @(x) ip.Results.G_fct; % If G_fct is a constant matrix, make it into a function
    end
    
    R_fct = ip.Results.R_fct;
    if isnumeric(R_fct)
        if isempty(R_fct)
            R_fct = @(x) ComputeG(ip.Results.r,x); % If R_fct is not supplied, use numerical differentiation
        else
            R_fct = @(x) ip.Results.R_fct; % If R_fct is a constant matrix, make it into a function
        end
    end
    
    assert(~isempty(res.sigma) || ~isempty(res.V), 'Either sigma or V must be supplied');
    assert(isempty(res.sigma) || length(res.sigma)==p, 'Dimension of sigma is incorrect');
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
        res.theta = theta_estim_fct(res.W);
    end
    res.theta = res.theta(:);
    res.h_theta = h(res.theta);
    k = length(res.theta); % Number of parameters
    
    res.G = G_fct(res.theta); % Moment function Jacobian at initial estimate
    assert(isequal(size(res.G),[p,k]), 'Jacobian of moment function has incorrect dimensions');
    res.R = R_fct(res.theta); % Transformation function Jacobian at initial estimate
    assert(size(res.R,2)==k, 'Jacobian of transformation function has incorrect dimensions');
    m = size(res.R,1); % Number of transformations of interest
    
    
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
                res.theta = theta_estim_fct(res.W);
            end
            
        else % Worst case analysis
            
            if m>1 % If there is more than one parameter of interest
                
                res.r_theta = nan(m,1);
                res.r_theta_se = nan(m,1);
                res.x_hat = nan(p,m);
                the_eye = eye(m);
                
                for im=1:m % Loop over each parameter of interest
                    % Recursive call to function with single parameter of interest
                    the_res = WorstCaseSE(h, res.mu, res.sigma, theta_estim_fct, ...
                                          'r', @(x) the_eye(im,:)*ip.Results.r(x), ...
                                          'W', res.W, 'theta', res.theta_init, ...
                                          'G_fct', res.G, 'R_fct', res.R(im,:), ...
                                          'opt', true, 'one_step', res.one_step, 'overid', false, ...
                                          'zero_thresh', ip.Results.zero_thresh);
                    res.r_theta(im) = the_res.r_theta;
                    res.r_theta_se(im) = the_res.r_theta_se;
                    res.x_hat(:,im) = the_res.x_hat;
                end
                
            else % Otherwise, compute worst-case optimum for the single parameter of interest
                
                % Worst-case optimal weighting and SE
                [res.x_hat, res.r_theta_se] = ComputeWorstCaseOpt_Single(res.sigma, res.G, res.R', ip.Results.zero_thresh);
                
                % Update point estimate
                if res.one_step % One-step update
                    res.r_theta = getOnestep(res.h_theta_init, [], res.x_hat, ip.Results.r(res.theta_init), res.mu);
                else % Full optimization update
                    % Weight matrix only puts weight on k moments with largest loadings in hat{x}
                    [~,I]=sort(abs(x_hat)); % Sort elements in hat{x} by magnitude
                    res.W(I(1:p-k),:) = 0;
                    res.W(:,I(1:p-k)) = 0;
                    res.theta = theta_estim_fct(res.W); % Update estimate
                end
                
            end
            
        end
        
        if isfield(res, 'theta') % Update Jacobians and moment estimate
            res.G = G_fct(res.theta);
            res.R = R_fct(res.theta);
            res.h_theta = h(res.theta);
        end
        
    end
    
    if ~isfield(res, 'r_theta')
        res.r_theta = ip.Results.r(res.theta); % Transformations of interest
    end
    
    
    %% SE
    
    if res.full_info % Full information analysis
        
        res.r_theta_varcov = ComputeCMDAvar(res.G,res.W,res.V,res.R); % Var-cov matrix
        res.r_theta_se = sqrt(diag(res.r_theta_varcov)); % SE
        
    else % Worst case analysis
        
        if ~res.opt % SE have already been computed above for worst-case optimal estimates
            
            res.x_hat = res.W*res.G*((res.G'*res.W*res.G)\res.R'); % Moment loadings
            res.r_theta_se = abs(res.x_hat)'*res.sigma; % SE
            
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
        the_res = WorstCaseSE(@(x) x, res.mu, res.sigma, [], ...
                              'r', @(x) M'*x, ...
                              'W', eye(p), 'V', res.V, 'theta', res.mu, ...
                              'G_fct', eye(p), 'R_fct', M', ...
                              'opt', false, 'overid', false, ...
                              'zero_thresh', ip.Results.zero_thresh); % Compute SE for M'*hat{mu}
        res.overid.se = the_res.r_theta_se;
        
        % Test statistics and p-values
        res.overid.t_stats = res.overid.errors./res.overid.se;
        res.overid.p_vals = 2*normcdf(-abs(res.overid.t_stats));
        
    end
    
end