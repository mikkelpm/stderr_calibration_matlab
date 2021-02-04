%% Compute worst case optimal W given V for different G matrices and parameters of interest
function [] = Test_ComputeWorstCaseOptW_Single()

  % Parameters of problem
  p    = 3; % Number of moments
  k    = 2; % Number of parameters

  % Set nonsingular V matrices for empirical moments

    % Construct different V matrices
    V = nan(p);
    V(find(eye(p))) = [9,4,1];


  % Element of G matrix that will vary to test how worst case optimal
  % changes as G matrix changes
  deltas  = -1:0.5:1;
  Ndeltas = length(deltas);

  % Preallocate to store output
  W      = nan(p,p,k,Ndeltas);
  stderr = nan(k,    Ndeltas);
  Vout   = nan(p,p,k,Ndeltas);
  x      = nan(p,k,  Ndeltas);
  z      = nan(p-k,k,Ndeltas);
  I      = eye(k);


  % Loop over G matrices and compute worst case optimal matrices
  for d = 1:Ndeltas
    delta = deltas(d);
    fprintf('Delta = %.2f\n', delta);

    % Set up G matrix
    G = [1      0;...
         delta  0;...
         0      1];

    % Compute worst-case optimal for each parameter of interest
    for n = 1:k
      [W(:,:,n,d), stderr(n,d), Vout(:,:,n,d), x(:,n,d), z(:,n,d)] = ComputeWorstCaseOptW_Single(V, G, I(:,n), 10e-8);
    end
  end


end
