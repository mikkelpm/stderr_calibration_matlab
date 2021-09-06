function loadings = get_moment_loadings(moment_jacob, weight_mat, transf_jacob)

    % Asymptotic loadings of minimum distance estimator on empirical moments
    
    loadings = weight_mat*moment_jacob*((moment_jacob'*weight_mat*moment_jacob)\(transf_jacob'));

end