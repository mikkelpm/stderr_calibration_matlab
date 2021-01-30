function [ bhat ] = qreg( y, x, tau )
%   bhat are the estimates
%   y is a vector of outcomes
%   x is a matrix with columns of explanatory variables
%   tau is a scalar for choosing the conditional quantile to be
%   estimated

n=size(x,1);
m=size(x,2);
% vectors and matrices for linprog
f=[tau*ones(n,1);(1-tau)*ones(n,1);zeros(m,1)];
Aeq=[eye(n),-eye(n),x];
beq=y;
lb=[zeros(n,1);zeros(n,1);-inf*ones(m,1)];
ub=inf*ones(m+2*n,1);

% Solve the linear programme
[bhat,~,~]=linprog(f,[],[],Aeq,beq,lb,ub);

% Pick out betas from (u,v,beta)-vector.
bhat=bhat(end-m+1:end);

end
