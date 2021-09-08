function param_estim = get_onestep(obj, moment_init, weight_mat, moment_jacob, param_init)

    % One-step estimation
    
    if isempty(weight_mat)
        subtr = moment_jacob'*(moment_init-obj.moment_estim);  
    else
        subtr = (moment_jacob'*weight_mat*moment_jacob)\(moment_jacob'*weight_mat*(moment_init-obj.moment_estim));
    end
    
    param_estim = param_init - subtr;
    
end

