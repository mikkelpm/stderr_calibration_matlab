% Do joint overid test based on checking all moments individually than
% seeing if any of the (Bonferroni adjusted) CIs don't contain zero
%
function [ret] = HTest_Overid(res, Nobs, p, alpha, wcopt)


  muhat = res.muhat;
  h     = res.h;
  G     = res.G;
  W     = res.W;
  V     = res.V;

  % The gap between the reduced-form moments and model-implied moments
  diff = muhat-h;

  % Calculate the variance of this that you will use for testing
  if wcopt
    % Worst case variance for that combination of moments
    Sigmas = abs(eye(p)-G*((G'*W*G)\(G'*W))) * sqrt(diag(V) / Nobs);
  else
    % Variance under the assumed (possibly incorrect or infeasible) V
    Sigmas  = sqrt( diag( (eye(p)-G*((G'*W*G)\(G'*W))) * V * (eye(p)-G*((G'*W*G)\(G'*W)))' ) / Nobs );
  end
  Sigmas0 = abs(Sigmas) < 1e-5;

  % Construct CIs for single hypotheses and joint overidentification
  cv_single   = abs(norminv(alpha/2));
  cv_joint    = abs(norminv(alpha/(2*p))); % Bonferroni adjusted
  ci_single_u = diff + cv_single*Sigmas;
  ci_single_d = diff - cv_single*Sigmas;
  ci_single   = [ci_single_d ci_single_u];
  ci_joint_u  = diff + cv_joint*Sigmas;
  ci_joint_d  = diff - cv_joint*Sigmas;
  ci_joint    = [ci_joint_d ci_joint_u];

  % Rejection of single hypotheses and joint hypothesis
  rej_single =      (0 > ci_single_u) + (0 < ci_single_d);
  rej_joint  = any( (0 > ci_joint_u)  + (0 < ci_joint_d) );

  % Correctly account for sigmas=0, i.e. when we can't test
  if wcopt && any(Sigmas0)
    rej_joint            = 0;
    rej_single(Sigmas0)  = 0;
    ci_single(Sigmas0,:) = NaN;
  end

  % Pack output into structure
  res.overid_diff        = diff;
  res.overid_stderr      = Sigmas;
  res.overid_cv_joint    = cv_joint;
  res.overid_ci_joint    = ci_joint;
  res.overid_rej_joint   = rej_joint;
  res.overid_cv_single   = cv_single;
  res.overid_ci_single   = ci_single;
  res.overid_rej_single  = rej_single;
  res.overid_nottestable = Sigmas0;

  ret = res;
end
