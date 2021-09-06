classdef MinDist
    
    % Minimum distance object used to compute worst-case standard errors and tests
    % Can also do standard full-information analysis
    
    % Reference:
    % Cocci, Matthew D. & Mikkel Plagborg-Moller, "Standard Errors for Calibrated Parameters"
    % https://scholar.princeton.edu/mikkelpm/calibration
    
    
    properties
        moment_fct
        moment_estim
        moment_num
        moment_fct_deriv
        moment_varcov
        full_info
        diag_only
        blockdiag_only
        moment_varcov_blocks
    end
    
    methods
        
        function obj = MinDist(moment_fct, moment_estim, varargin)
            
            % Class constructor
            
            if nargin>0
                
                % Parse inputs
                ip = inputParser;
                addRequired(ip, 'moment_fct', @(x) isa(x, 'function_handle'));
                addRequired(ip, 'moment_estim', @isnumeric);
                addParameter(ip, 'moment_se', [], @isnumeric);
                addParameter(ip, 'moment_varcov', [], @isnumeric);
                addParameter(ip, 'moment_fct_deriv', [], @(x) isa(x, 'function_handle') | isempty(x));
                parse(ip, moment_fct, moment_estim, varargin{:}); 
                
                % Store values
                obj.moment_fct = ip.Results.moment_fct;
                obj.moment_estim = ip.Results.moment_estim(:);
                obj.moment_num = length(ip.Results.moment_estim);
                obj.moment_fct_deriv = deriv(ip.Results.moment_fct_deriv, ip.Results.moment_fct);
                
                % Var-cov matrix
                assert(isempty(ip.Results.moment_se) | isempty(ip.Results.moment_varcov), ...
                       'Either "moment_se" or "moment_varcov" must be supplied, but not both');
                if ~isempty(ip.Results.moment_varcov)
                    obj.moment_varcov = ip.Results.moment_varcov;
                else
                    obj.moment_varcov = nan(obj.moment_num);
                    obj.moment_varcov(eye(obj.moment_num)==1) = ip.Results.moment_se.^2; % Diagonal known
                end
                
                % Check inputs
                assert(size(obj.moment_varcov,1)==obj.moment_num & size(obj.moment_varcov,2)==obj.moment_num, ...
                       'Dimension of "moment_se" or "moment_varcov" is wrong');
                assert(all(~isnan(obj.moment_estim)), ...
                       'Wrong input type for "moment_estim"');
                assert(all(~isnan(diag(obj.moment_varcov))), ...
                       'Wrong input type for "moment_se" or "moment_varcov"');
                assert(all(diag(obj.moment_varcov)>=0), ...
                       'SE for each individual moment must be nonnegative');
                assert(all(obj.moment_varcov==obj.moment_varcov' | isnan(obj.moment_varcov),'all'), ...
                       '"moment_varcov" must be symmetric');
                assert(obj.moment_num>=1, ...
                       '"moment_estim" must have at least one element');
               
               % Determine type of var-cov matrix
               obj.full_info = ~any(isnan(obj.moment_varcov),'all'); % Full info
               obj.diag_only = (~obj.full_info & all(isnan(obj.moment_varcov(eye(obj.moment_num)==0)),'all')); % Only diagonal is known
               
               % Check if var-cov has known block diagonal (and unknown everywhere else)
               obj.blockdiag_only = false;
               if ~obj.full_info && ~obj.diag_only
                   i = 1;
                   block_inds = [];
                   while i<=obj.moment_num % Loop through purported blocks
                       the_block = i-1 + find(~isnan(obj.moment_varcov(i,i:end)));
                       if ~all(diff(the_block)==1)
                           return; % Can't be block diagonal
                       end
                       block_inds = [block_inds {the_block}];
                       i = max(the_block)+1;
                   end
               
                   % Check that the block diagonal indeed has non-NaN values, with NaN's outside blocks
                   block_ones = cell(1,length(block_inds));
                   for j=1:length(block_ones)
                       block_ones{j} = ones(length(block_inds{j}));
                   end
                   block_bool = blkdiag(block_ones{:});
                   obj.blockdiag_only = (all(~isnan(obj.moment_varcov(block_bool==1))) & all(isnan(obj.moment_varcov(block_bool==0))));
                   if ~obj.blockdiag_only
                       return;
                   end
                   
                   % If block diagonal, extract useful information for later
                   obj.moment_varcov_blocks = struct;
                   obj.moment_varcov_blocks.inds = block_inds; % Indices of blocks
                   obj.moment_varcov_blocks.num = length(block_inds);
                   obj.moment_varcov_blocks.varcov = cell(1,length(block_inds));
                   obj.moment_varcov_blocks.chol = cell(1,length(block_inds));
                   for j=1:length(block_inds)
                       obj.moment_varcov_blocks.varcov{j} = obj.moment_varcov(block_inds{j},block_inds{j}); % Var-cov for each block
                       obj.moment_varcov_blocks.chol{j} = chol(obj.moment_varcov_blocks.varcov{j},'lower'); % Cholesky factors
                   end
               
               end
                
            end
            
        end
        
    end
end
    
    