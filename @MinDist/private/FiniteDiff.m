% For function h(theta): Rk -> Rp, compute G = dh/dtheta' by finite diff
function [G] = FiniteDiff(h, theta)

  h_theta   = h(theta);
  p         = length(h_theta);
  k         = length(theta);
  delta     = 0.01;   % Size of param perturbation. Will shrink until computed deriv converges
  thresh    = 1e-6;   % Threshold for convergence
  converged = 0;      % Whether computed derivative has converged
  iter      = 0;      % Flag the first run trhough
  maxabspctchg_old = Inf;
  while ~converged

    % Finite difference perturbation amount
    Delta = delta*eye(k);

    % ith column of theta_fwd is theta, but with its ith element
    % perturbed forward by delta. theta_bwd similar.
    theta_fwd = repmat(theta, 1, k) + Delta;
    theta_bwd = repmat(theta, 1, k) - Delta;

    % Compute h under each perturbed theta value
    h_fwd = nan(p, k);
    h_bwd = nan(p, k);
    for i_ = 1:k
      h_fwd(:,i_) = h(theta_fwd(:,i_));
      h_bwd(:,i_) = h(theta_bwd(:,i_));
    end
    h_ = repmat(h_theta, 1, k);

    % Compute estimate of G
    G_fwd = (h_fwd - h_)/delta;
    G_bwd = (h_ - h_bwd)/delta;
    G_ctr = (h_fwd - h_bwd)/(2*delta);
    G_new = G_ctr;

    % Check for convergence in derivative if not iter=0
    if iter
      abspctchg    = abs((G_new - G_old) ./ G_old);
      maxabspctchg = max(abspctchg(:));
      wrongdir     = (maxabspctchg > maxabspctchg_old);
      converged    = (maxabspctchg < thresh) || wrongdir;

      if ~converged
        if isnan(max(abspctchg(:)))
          error('Error: max abs pct chg = NaN')
        end
%         fprintf('\nNorm = %f\n', max(abspctchg(:)))
%         fprintf('\nShrinking to delta = %f\n', 0.1*delta)
      end
      maxabspctchg_old = maxabspctchg;
    end

    % Update delta, G_old, iter in case needed for next iteration
    delta            = 0.1*delta;
    G_old            = G_new;
    iter             = iter+1;
  end
  if wrongdir
    G = G_old;
  else
    G = G_new;
  end

end
