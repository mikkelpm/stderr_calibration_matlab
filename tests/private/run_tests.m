function [res, res_eff] = run_tests(obj, G, mu, theta, tol)

    % Run a bunch of simple tests on output

    % Estimation output
    res = obj.fit('opt_init', zeros(size(theta)), 'eff', false);
    assert(isequal(size(res.moment_loadings), size(G)));
    assert(norm(res.moment_fit - mu)<tol);
    assert(norm(res.moment_jacob - G)<tol);
    assert(norm(res.estim - theta)<tol);
    assert(all(res.estim_se>=0));
    
    % Test output
    test_res = obj.test(res);
    assert(test_res.joint_stat>=0);
    overid_res = obj.overid(res);
    assert(norm(overid_res.moment_error)<tol);
    assert(overid_res.joint_stat>=0);
    
    % Efficient estimation
    res_eff = obj.fit('opt_init', zeros(size(theta)));
    assert(norm(res_eff.estim - theta)<tol);

end