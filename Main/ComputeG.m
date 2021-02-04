function [G] = ComputeG(h, theta)
  G = FiniteDiff(h, theta);
end
