%% Setup

  % Make sure cvx package loads upon startup. See Readme for details

  % Directory with test examples
  addpath Test/

  % Directory with code to solve the semidefinite programming problems and
  % worst-case variance problems
  addpath Main/

  % Directory with supporting code to help out the functions in Main
  addpath Supporting/


%% Testing

% Test the program that solves the SDP program max x'Vx for V incomplete
Test_SDP(0)   % Do not do speed tests
%Test_SDP(1)  % Do speed test. Takes a long time.


% Test the program that solves for the worst-case optimal W matrix when
% testing a hypothesis about a single parameter using the quantile
% regression approach approved in the paper. Do this for multiple G
% matrices.
Test_ComputeWorstCaseOptimal_Single
