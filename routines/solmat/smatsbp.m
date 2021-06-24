% smatsbp.m 
%
% [Q, G, unique] = smatsbp(A0,A1,B0,D0)
%
% Sytem is given as 
%
%  A0 y(t) = A1 y(t-1) + B0 E_t y(t+1) + D0 w(t) + D2 E_t w(t+1)
%
%     w(t) = Gamma w(t-1) + e(t)
%
%  and gives solution of the form
%  
%  y(t) = Q y(t-1) + G w(t)
%
% Dimensions: y(t) is nvars x 1
%             w(t) is nshocks x 1
%             e(t) is nshocks x 1
%
% Based on code by Binder and Pesaran (1995,1997)
%
% Date created: 29/12/2010
% Modified: Mariano Kulish
% [add notes here]


function [Q, G, unique, rc, iter] = smatsbp(A0,A1,B0,D0,D2,Gamma);
%%
N = 300; % Initial Forecasting Horizon (See Binder and Pesaran, 1995)
unique = 1; % equal to 1 if unique solution exists and 0 otherwise.

A = A0\A1;
B = A0\B0;


[dim1,dim2] = size(A);
[dim3,dim4] = size(Gamma);
% 
% Compute Matrix Q Using Brute-Force Iterative Procedure
Q = eye(dim1);          % Initial Condition
F = eye(dim1);          % Initial Condition
eps1 = 10^(-6);         % Convergence Criterion for F
eps2 = 10^(-6);         % Convergence Criterion for Q
crit1 = 1; crit2 = 1;   % Initial Conditions
iter = 0;
while crit1 >= eps1 || crit2 >= eps2
   Fi = (eye(dim1)-B*Q)\B;
   Ci = (eye(dim1)-B*Q)\A;
   rc  = rcond(eye(dim1)-B*Q);
   crit1 = max(max(abs(Fi-F))); crit2 = max(max(abs(Ci-Q)));
   Q = Ci; F = Fi;
   iter = iter+1;
   
   % Returns (true) if either input, or both, evaluate to true, and (false) if they do not.
   if iter > 500 || rc < 1.0e-15, 
%       disp(' Did not converge after 500 iterations')
%       disp(' or (I-BQ) badly conditioned.')
%       disp(' INDETERMINACY OR NON-EXISTENCE')
      unique = 0;
      Q = []; G = []; return
   end
end

%% Use Recursive Method of Binder and Pesaran (1995) to Compute the 
% Forward Part of the Solution - Determine N
eps3 = 10^(-6);        % Convergence Criterion
i = 0;
aux3a = zeros(dim1,dim3);
aux3b = zeros(dim1,dim3);
while i <= N
   fp1 = F^i/(eye(dim1)-B*Q)/A0*D0*Gamma^i;
   fp2 = F^i/(eye(dim1)-B*Q)/A0*D2*Gamma^(i+1);
   aux3a = fp1+aux3a;
   aux3b = fp2+aux3b;
   i = i+1;
end
crit3 = max(max(abs(fp1+fp2)));

while crit3 > eps3
   N = N+1;
   fp1 = F^N/(eye(dim1)-B*Q)/A0*D0*Gamma^N;
   fp2 = F^N/(eye(dim1)-B*Q)/A0*D2*Gamma^(N+1);
   aux3a = fp1+aux3a;
   aux3b = fp2+aux3b;
   crit3 = max(max(abs(fp1+fp2)));
end
G = aux3a+aux3b;


