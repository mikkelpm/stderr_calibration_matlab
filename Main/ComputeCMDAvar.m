% Formula for asymptotic variance in classical minimum distance
function [Omega] = ComputeCMDAvar(G,W,V,lambda)
  Omega = lambda'*((G'*W*G)\G')*W*V*W*G*((G'*W*G)\lambda);
end



