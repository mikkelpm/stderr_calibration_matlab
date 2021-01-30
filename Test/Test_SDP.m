% Testing the time to solve a semidefinite programming problem using CVX
function [] = Test_SDP(speedtest)


% Simple Example 1
Vtilde      = nan(3,3);
Vtilde(1,1) = 1;
Vtilde(2,2) = 1;
Vtilde(3,3) = 1;
lambda = ones(3,1);
A      = lambda*lambda';
[Vtilde, cvx_optval, cvx_status] = SolveVarianceSDP(A, Vtilde)


% Simple Example 2
Vtilde      = nan(3,3);
Vtilde(1,1) = 1;
Vtilde(2,2) = 2;
Vtilde(3,3) = 3;
lambda = ones(3,1);
A      = lambda*lambda';
[Vtilde, cvx_optval, cvx_status] = SolveVarianceSDP(A, Vtilde)


% Simple Example 3
Vtilde      = nan(3,3);
Vtilde(1,1) = 1;
Vtilde(2,2) = 2;
Vtilde(3,3) = 3;
Vtilde(1,2) = 1;
lambda = ones(3,1);
A      = lambda*lambda';
[Vtilde, cvx_optval, cvx_status] = SolveVarianceSDP(A, Vtilde)


%% Check computation times
%
% Maximize x'Cx s.t. only diagonals of C known (and equal to 1) for
% randomly generated x. Repeat this many times.
%
% Do this for different sizes of the C and x. Store and plot results.
%
outerprod = @(X) X*X';
if speedtest
  %profile on

  % Number of parameters
  ks = 1:20;
  Nk = length(ks);

  % Number of moments to loop over. Always assume p > k
  pmax = 40;
  Np   = pmax;

  % Number of reps for each (k,p) combination
  Nrep = 100;

  % For storing computation times
  time = nan(Nk, Np, Nrep);

  % Loop and compute
  for k = ks
    fprintf('\nNumber of Parameters: %d\n', k)

    for p = (k+1):pmax
      fprintf('\tNumber of Moments: %d\n', p)
      C = nan(p,p);
      C( find(diag(ones(p,1))) ) = 1;
      for n = 1:Nrep
        tic();
        SolveVarianceSDP(outerprod(randn(p,k)), C);
        time(k,p,n) = toc();
      end
    end
  end
  save ComputationTimes_new.mat
  %profile off


  % Plot
  load ComputationTimes.mat
  figure()
  for p = ps
    for k = p:30
      subplot(Nk, Np, (k-1)*Np + p)
      histogram(squeeze(time(k,p,:)))
      title(sprintf('(k,p)=(%d,%d)', k, p))
    end
  end

end




end
