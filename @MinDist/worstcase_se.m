function [se, varcov] = worstcase_se(obj, moment_loadings)

    % Worst-case standard errors and corresponding var-cov matrix
    % for linear combination of moments
    
    
    estim_num = size(moment_loadings,2);
    if estim_num>1 % If more than one parameter of interest, handle them separately
        se = nan(estim_num,1);
        varcov = cell(1,estim_num);
        for i=1:estim_num
            [se(i),varcov{i}] = obj.worstcase_se(moment_loadings(:,i));
        end
        return;
    end
    
    if obj.diag_only % Only diagonal is known
        
        moment_se = sqrt(diag(obj.moment_varcov));
        se = abs(moment_loadings)'*moment_se;
        aux = sign(moment_loadings).*moment_se;
        varcov = aux*aux';
        
    elseif obj.blockdiag_only % Only block diagonal is known
        
        nb = obj.moment_varcov_blocks.num; % Number of blocks
        var_blocks = nan(nb,1);
        varcov_left = nan(nb,1);
        varcov_block = cell(1,nb);
        
        for j=1:nb
            the_inds = obj.moment_varcov_blocks.inds{j};
            the_loadings_block = moment_loadings(the_inds,:);
            the_chol = obj.moment_varcov_blocks.chol{j};
            aux = obj.moment_varcov_blocks.varcov{j} * the_loadings_block;
            var_blocks(j) = max(the_loadings_block' * aux, 1e-10); % Avoid exact zeros (when loadings are zero)
            varcov_left(the_inds) = aux/sqrt(var_blocks(j));
            varcov_block{j} = the_chol - aux*(the_loadings_block'*the_chol)/var_blocks(j);
        end
        
        se = sum(sqrt(var_blocks));
        varcov_aux = [varcov_left, blkdiag(varcov_block{:})];
        varcov = varcov_aux*varcov_aux';
        
    else % General knowledge of var-cov matrix
        
        % Solve semidefinite programming problem
        [vari, varcov] = obj.solve_sdp(moment_loadings*moment_loadings');
        se = sqrt(vari);
        
    end

end