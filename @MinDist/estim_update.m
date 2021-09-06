function [estim, transf_jacob, moment_fit, moment_jacob] = ...
        estim_update(obj, param_estim, transf, transf_deriv, ...
                     estim, transf_jacob, moment_fit, moment_jacob)
     
    % Update estimated parameter transformation and moments, including their Jacobians
    % Avoids recomputing quantities if they're already supplied
    
    if isempty(estim)
        estim = transf(param_estim);
    end
    
    if isempty(transf_jacob)
        transf_jacob = transf_deriv(param_estim);
    end
    
    if isempty(moment_fit)
        moment_fit = obj.moment_fct(param_estim);
    end
    
    if isempty(moment_jacob)
        moment_jacob = obj.moment_fct_deriv(param_estim);
    end

end