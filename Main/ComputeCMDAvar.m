% Formula for asymptotic variance in classical minimum distance
function [Omega] = ComputeCMDAvar(G,W,V,M)
  Omega = M*((G'*W*G)\G')*W*V*W*G*((G'*W*G)\M');
end



