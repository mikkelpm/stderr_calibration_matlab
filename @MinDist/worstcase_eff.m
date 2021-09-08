function [se, moment_loadings, weight_mat_new] = worstcase_eff(obj, moment_jacob, transf_jacob, weight_mat)

    % Compute worst-case efficient moment loadings and weight matrix
    % See main paper for explanation


    % Set up median regression as described in paper
    [p,k]  = size(moment_jacob);
    GpG = moment_jacob'*moment_jacob;
    Y = moment_jacob*(GpG\(transf_jacob'));
    moment_jacob_perp = null(moment_jacob');
    X = -moment_jacob_perp;

    if obj.diag_only % Only diagonal is known
        
        % Run median regression
        moment_se = sqrt(diag(obj.moment_varcov));
        z = qreg(moment_se.*Y, moment_se.*X, 0.5);
        moment_loadings = Y-X*z;
        se = moment_se'*abs(moment_loadings);

        % Weight matrix puts weight on only k moments
        [~,sort_inds] = sort(abs(moment_loadings));
        if ~isempty(weight_mat)
            weight_mat_new = weight_mat;
        else
            weight_mat_new = eye(p);
        end
        weight_mat_new(sort_inds(p-k+1:end),:) = 0;
        weight_mat_new(:,sort_inds(p-k+1:end)) = 0;
        
    else % General case
        
        % Solve nested optimization
        opts = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'Display', 'notify');
        [z_opt, vari] = fminunc(@(z) objective(z,Y,X,obj), zeros(p-k,1), opts);
        moment_loadings = Y-X*z_opt;
        se = sqrt(vari);
        
        % Weight matrix
        aux1 = (transf_jacob'*z_opt')/(transf_jacob*(GpG\(transf_jacob')));
        W = @(delta) [eye(k), aux1; aux1', delta*eye(p-k)];
        % Determine delta such that W(delta) is positive definite
        opts2 = optimoptions('fsolve', 'Display', 'none');
        delta_pd = fsolve(@(delta) eigs(W(delta),1,'smallestabs')-0.01, 0, opts2);
        aux2 = [moment_jacob, moment_jacob_perp];
        weight_mat_new = aux2 * W(delta_pd) * aux2';
        
    end
  
end


function [f,grad] = objective(z, Y, X, obj)
    
    % Objective function and gradient for general case
    
    resid = Y-X*z;
    [se, varcov] = obj.worstcase_se(resid);
    f = se*se;
    grad = -2*(X'*varcov*resid);

end