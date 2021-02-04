% Computing the one-step estimator for classical minimum distance
%
% h0      h() at initial consistent estimate, theta0
%
% W       Weighting matrix
%
% G       Initial consistent estimate of G
%
% theta0  Initial constistent estimate of theta
%
% muhat   Empirical moments
%
function [theta] = getOnestep(h0, W, G, theta0, muhat)
  subtr = (G'*W*G)\(G'*W*(h0-muhat));
  theta = theta0 - subtr;
end

