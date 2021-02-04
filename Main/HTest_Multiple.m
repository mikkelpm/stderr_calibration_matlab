% Under the null
%
% sqrt(n)(R thethat - r) -> N(0, R (G'WG)^{-1} G'W V WG (G'WG)^{-1} R')
%
% Let H = R (G'WG)^{-1} G'W
%
function [reject, cv, pvalue] = HTest_Multiple(n, thetahat, R, r, G, W, V, S)

  % Compute test stat
  stat = n*(R*thetahat-r)'*(S\(R*thetahat-r));

  % Get the Chi2_1 critical value, which is part of the full cv
  cv_chi2 = chi2inv(0.95, 1);

  % Find the number that we have to blow up the critical value by
  H = R*((G'*W*G)\(G'*W));
  inflate = SolveVarianceSDP((H'*(S\H))', V);

  % Determine critical value
  cv = inflate*cv_chi2;

  % Determine accept/reject
  reject = (stat > cv);

  % Determine p-value if test stat bigger than cv (in which case
  % computing a p-value is valid, otherwise it's not)
  pvalue = NaN;
  if reject
    pvalue = chi2cdf(stat/inflate, 1, 'upper');
  end

end
