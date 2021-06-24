% SMATS.m
%
%[S0, S1, S2, S3, Sy, Sf, Sz, unique] = Smats(C, GAM0, GAM1, PSI, PPI);
%
% System given as
%        GAM0*y(t)=GAM1*y(t-1)+C+PSI*e(t)+PPI*eta(t),
% with e(t) an exogenous variable process and eta(t) being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=S0 + S1*y(t-1)+S2*e(t)+S3*sunspot(t)+ Sy sum(j=1..inf, Sf^(j-1) Sz e(t+j))  .

% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.

% Sims (2002)
%
% Edited by Callum Jones, jonescallum@gmail.com
%-------------------------------------------------------------------------%

function [S0, S1, S2, S3, Sy, Sf, Sz, unique] = smats(C, GAM0, GAM1, PSI, PPI);


%% Start Calculations

eu=[0;0];

realsmall=1e-6;

[n k] = size(PPI);
[n l] = size(PSI);

[A B Q Z ~]=qz(GAM0,GAM1);
[A B Q Z ~]=ordqz(A,B,Q,Z,'udo'); % Added to original code and commented out lines 56 to 59 below

div = 1.01;
nunstab=0;
zxz=0;

for i=1:n % Determines the dimension of the unstable block of the model. The size of W_1(t)is n-m 
       if abs(A(i,i)) > 0
          divhat=abs(B(i,i))/abs(A(i,i));
          if 1+realsmall<divhat && divhat<=div
             div=.5*(1+divhat);
          end
       end
   nunstab = nunstab + (abs(B(i,i))>div*abs(A(i,i)));
   if abs(A(i,i))<realsmall && abs(B(i,i))<realsmall % Checking coincident zeros
      zxz=1;
   end
end

div ;
nunstab;
m = nunstab;

% if ~zxz % If no coincident zeros do QZDIV.m
%    [A B Q Z]= qzdiv(div,A,B,Q,Z);
% end
% vin = 0;
%keyboard
if m>k 
    disp('No stable solution exists') ;
    unique = -1 ;
    S0 = []; S1=[]; S2=[]; S3=[]; Sy=[]; Sf=[]; Sz=[];
    return
end

%keyboard

%%
gev=[diag(A) diag(B)]; % The ratios diag(B)/diag(A) are the generalized eigenvalues of GAM1 and GAM0

if zxz
   disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
   eu=[-2;-2];
   % correction added 7/29/2003.  Otherwise the failure to set output
   % arguments leads to an error message and no output (including eu).
   %disp('Generalized eigenvalues')
   %disp(sort(abs(gev(:,2) ./ gev(:,1))));
   %disp('Dimension of eta is:')
   %disp(k)
   S0 = []; S1=[]; S2=[]; S3=[]; unique=0; Sy=[]; Sf=[]; Sz=[];
   return
end
%%

Q1=Q(1:n-nunstab,:);
Q2=Q(n-nunstab+1:n,:);
Z1 = Z(:,1:n-nunstab)';
Z2 = Z(:,n-nunstab+1:n)';
A22 = A(n-nunstab+1:n,n-nunstab+1:n); % Lambda22 in LS (2003)
B22 = B(n-nunstab+1:n,n-nunstab+1:n); % Omega22 in LS (2003)
etawt=Q2*PPI;

% Input Matrices for Indeterminate Solution 

A11 = A(1:n-nunstab,1:n-nunstab);
A12 = A(1:n-nunstab,n-nunstab+1:n);


% SVD of ETAWT 
[U,D,V] = svd(etawt);
md = min(size(D));
bigev = find(diag(D(1:md,1:md))>realsmall); % Partitions the SVD according to the number of non-zero singular values (LS, JEDC - Eq9, p277)
U1 = U(:,bigev); %LS: U.1 ;
V1 = V(:,bigev); %LS: V.1
D11 = D(bigev,bigev); %LS: D11
r = size(D11,1);
V2 = V(:,r+1:end);  
U2 = U(:,r+1:end);  

eu(1) = length(bigev)>=nunstab;
%%
[Uhat,Dhat,Vhat] = svd(Q1*PPI);
md = min(size(Dhat));
bigev = find(diag(Dhat(1:md,1:md))>realsmall);
Uhat1 = Uhat(:,bigev);
Vhat1 = Vhat(:,bigev);
Dhat11 = Dhat(bigev,bigev);
rr = size(Dhat11,1);
Uhat2 = Uhat(:,rr+1:end); 
Vhat2 = Vhat(:,rr+1:end);
%%

estimate = 1; % Set to 0 when not estimating. 

if eu(1)
    if estimate == 0
    disp('Solution Exists')
    end
else
   disp('Solution Does Not Exist')
   eu ;
   %disp('Generalized eigenvalues')
   %disp(sort(abs(gev(:,2) ./ gev(:,1))));
   %disp('Dimension of eta is:')
   %disp(k)
   S0 = []; S1=[]; S2=[]; S3=[]; unique=0; Sy=[]; Sf=[]; Sz=[];
   return
end

% if  isempty(V2)
% 	unique=1;
%     disp('The solution is unique');
%         
% else
%     unique=0;
% 	disp('The solution is indeterminate');
% end
% 
% if unique
%    %disp('solution unique');
%    eu(2)=1;
% else
%    fprintf(1,'Indeterminacy.  %d loose endog errors.\n',k-r);
%    disp('Generalized eigenvalues')
%    disp(gev);
%    disp(sort(abs(gev(:,2) ./ gev(:,1))));
% 
% end

%% Check for Uniqueness
if  isempty(Vhat1)
	unique=1;
else
	unique=norm(Vhat1-V1*V1'*Vhat1)<realsmall*n;
end

if unique
    if estimate == 0
       disp('Solution Unique');
    end
   %disp('Generalized eigenvalues')
   %disp(sort(abs(gev(:,2) ./ gev(:,1))));
   eu(2)=1;
else
   fprintf(1,'Indeterminacy.  %d loose endog errors.\n',size(Vhat1,2)-size(V1,2));
   %disp('Generalized eigenvalues')
   %disp(sort(abs(gev(:,2) ./ gev(:,1))));
   %disp(gev);
   %md=abs(diag(a))>realsmall;
   %ev=diag(md.*diag(a)+(1-md).*diag(b))\ev;
   %disp(ev)
%   return;
end

%%
PHI = (U1*(D11\V1')*Vhat1*Dhat11*Uhat1')';

tmat = [eye(n-nunstab) -PHI];
% Calculate matrices which would appear in the solution
G0 = [ tmat * A ; 
       zeros(m,n-m) eye(m)];
   
G1 = [ tmat * B ; 
       zeros(m,n)]; 

G2 = [zeros(n-m,m);-eye(m)];
   
G0I=inv(G0);

H = Z*G0I;

%%

if ~unique
    p = 1;  % Set dimension of the sunspot vector
    m_1 = 0;%input('Enter 0 for M1 of zeros and 1 for least squares M1: ');
    m_2 = 1;%input('Enter 0 for M2 of zeros and 1 for M2 of ones: ');
    
    if m_1
       K1 = Q1*PPI-PHI*Q2*PPI;
       K2 = -V1*inv(D11)*U1'*Q2*PSI;
       AA = K1*K2;
       BB = K1*V2; 
       M1 = inv(BB'*BB)*BB'*(-AA);
    else
       M1 = zeros(k-r,l);
    end
    if m_2
       M2 = ones(k-r,p);
    else
       M2 = zeros(k-r,p);
    end
end

%% Calculate Solution Matrices
S0 = real(H*[Q1 - PHI * Q2 ; inv(A22-B22)*Q2 ]*C);
S1 = real(H*G1*Z');
Sy = real(H*G2);
Sf = real(B22\A22);
Sz = real(B22\(Q2*PSI));

if unique  
   S2 = H * [ Q1 - PHI * Q2 ; zeros(m,n) ] * PSI;
   S2 = real(S2);
   
   S3 = zeros(0);
else
   S2 = H*[Q1-PHI*Q2 ; zeros(m,n)]*PSI + H*[(Q1*PPI-PHI*Q2*PPI)*(-V1*inv(D11)*U1'*Q2*PSI+V2*M1); zeros(m,l)];
   S2 = real(S2);
   
   S3 = H * [ (Q1*PPI - PHI*Q2*PPI)*V2*M2 ; zeros(m,p)];
   S3 = real(S3);
end




