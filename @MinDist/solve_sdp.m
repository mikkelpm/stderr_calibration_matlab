function [cvx_optval, Vtilde, cvx_status] = solve_sdp(obj, A)

    % Optimize over symmetric positive semidefinite matrices, solving
    %
    %   max_Vtilde  trace(A * Vtilde)
    %
    %   s.t. Vtilde matches V (the moment var-cov matrix) wherever V ~= NaN
    %
    % Requires cvx package to be loaded to solve this semidefinite
    % programming problem.

    p = obj.moment_num;
    V = obj.moment_varcov;

    cvx_begin quiet
      % Declare Vtilde, symmetric psd
      variable Vtilde(p,p) symmetric;
      Vtilde == semidefinite(p);

      % Set the entries of Vtilde specified as nonmissing in V
      Vtilde(~isnan(V)) == V(~isnan(V));

      % Solve
      maximize( trace(A*Vtilde) )
    cvx_end

end
