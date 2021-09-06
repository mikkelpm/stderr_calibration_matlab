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
  if isempty(W)
    subtr = G'*(h0-muhat);  
  else
    subtr = (G'*W*G)\(G'*W*(h0-muhat));
  end
  theta = theta0 - subtr;
end

