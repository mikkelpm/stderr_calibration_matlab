% Compute worst-case optimal weighting matrix for testing hypotheses
% about a single parameter in classical minimum distance problem, using
% quantile regression approach outlined in paper.
%
% V               Var matrix with NaNs off diagonal reflectin unknown elements.
%
% G               Consistent estimator for G matrix
%
% lambda          lambda giving the linear combination lambda*theta that
%                 we want to test a hypothesis about.
%
% zero_thresh     Ensuring that Gperp grabs the right eigenvectors, i.e.
%                 does not grab eigenvectors with nonzero-but-near-zero
%                 eigenvalues (arising from numerical precision issues)
%
function [x_hat, stderr, z] = ComputeWorstCaseOpt_Single(sigma, G, lambda, zero_thresh)


  % Set up median regression as described in paper
  p      = length(sigma);
  Y      = diag(sigma)*G*((G'*G)\lambda);
  [Q, D] = eig(eye(p)-G*((G'*G)\G'));
  Gperp  = Q(:, abs(diag(D)) > zero_thresh);
  X      = -diag(sigma)*Gperp;

  % Run median regression
  z = qreg(Y, X, 0.5);
  res = Y-X*z;

  % Standard errors
  stderr = sum(abs(res));
  
  % Moment loadings
  x_hat = res./sigma;
  
end
