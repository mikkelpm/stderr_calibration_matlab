% Optimize over symmetric positive semidefinite matrices, solving
%
%   max_Vtilde  trace(A' Vtilde)
%
%   s.t. Vtilde matches V (the input argument) wherever V ~= NaN
%
% Requires cvx package to be loaded to solve this semidefinite
% programming problem.
%
% Inputs
% ------
% V   (p x p) variance matrix. Use NaN for unknown elements (that this
%     problem will solve for and fill in). Filled-in values in V will
%     be used as constraints, so Vtilde(m,n)=V(m,n) by construction if
%     V(m,n) =/= NaN
%
% A   Weights on the elements of V
%
function [Vtilde, cvx_optval, cvx_status] = SolveVarianceSDP(A, V)

p = length(V);

cvx_begin quiet
  % Declare Vtilde, symmetric psd
  variable Vtilde(p,p) symmetric;
  Vtilde == semidefinite(p);

  % Set the entries of Vtilde specified as nonmissing in V
  Vtilde(~isnan(V)) == V(~isnan(V));

  % Solve
  maximize( trace(A'*Vtilde) )
cvx_end


end
