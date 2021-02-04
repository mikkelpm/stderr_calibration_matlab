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
function [W, stderr, Vout, x, z] = ComputeWorstCaseOptimal_Single(V, G, lambda, zero_thresh)

  % Differs from Matlab's sign fcn in that this has sign(0) = 1, rather
  % than sign(0)=0.
  %
  % Need the former to correctly implement the formulas from the paper
  % for the worst case optimal.
  mysign = @(x) (x >= 0) + -1*(x < 0);

  % Set up median regression as described in paper
  p      = length(V);
  Y      = diag(sqrt(diag(V)))*G*((G'*G)\lambda);
  [Q, D] = eig(eye(p)-G*((G'*G)\G'));
  Gperp  = Q(:, abs(diag(D)) > zero_thresh);
  X      = -diag(sqrt(diag(V)))*Gperp;

  % Run median ressression
  z = qreg(Y, X, 0.5);



  % Construct optimal weight matrix using median regression output
  %
  % If the eigenvalue of W is negative (hence W not psd) from numerical
  % precision issues in quantile regression, can include additional term
  % in expression that will not change the worst case optimal x(W) and
  % parameter estimate, but will ensure that W is positive semidefinite.
  %
  W = (G*G') + (Gperp*z*lambda'*G' + G*lambda*z'*Gperp') / (lambda'*((G'*G)\lambda));

  delta = 0.000001;
  while min(eig(W)) < 0
    delta = 10*delta;
    W = (G*G') + delta*Gperp*Gperp' + (Gperp*z*lambda'*G' + G*lambda*z'*Gperp') / (lambda'*((G'*G)\lambda));
  end


  % Construct V
  x = W*G*inv(G'*W*G)*lambda;
  s = sqrt(diag(V)).*mysign(x);
  Vout = s*s';

  % Compute standard errors
  stderr = sum(abs(Y-X*z));
end
