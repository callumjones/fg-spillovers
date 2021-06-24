function out = bp2smats(M_,A0,A1,B0,D0)
%
% Goes from Binder-Peseran matrices to Sims matrices.
%
% Solves using smats.
% 
% Callum Jones
% jonescallum@gmail.com

%%

M_ind  = M_.lead_lag_incidence ;
n_     = M_.endo_nbr+1 ;       
l_     = M_.exo_nbr ;        
n2     = sum(M_ind(3,:)>0) ; 
n1     = n_-n2-1 ;             
k_     = n2 ;

bp2smats_ord = 1:n_ ;
bp2smats_ord(end) = [] ;
bp2smats_ord = [ ...
    bp2smats_ord(M_ind(3,:)==0) ...
    bp2smats_ord(M_ind(3,:)>0)] ;

[~,smats2bp_ord] = sort(bp2smats_ord) ;

%ord_nonE = SET.variable.bp2smats_ord(1:n1) ;
ord_E = bp2smats_ord(n1+1:n1+n2) ;

%%

GAM0t = zeros(n1+n2,n1+n2+k_);
GAM1t = zeros(n1+n2,n1+n2+k_);
PSIt  = zeros(n1+n2,l_);
Ct    = zeros(n1+n2,1);

A0_ = A0(1:n_-1,bp2smats_ord) ;
B0_ = B0(1:n_-1,ord_E) ;
A1_ = A1(1:n_-1,bp2smats_ord) ;

GAM0t(1:n_-1,:) = [A0_ -B0_] ;
GAM1t(1:n_-1,:) = [A1_ zeros(n_-1,n2)] ;
Ct(1:n_-1,1)    = A0(1:n_-1,end) ;
PSIt(1:n_-1,:)  = D0(1:n_-1,:) ;

GAM0   = [GAM0t; [zeros(k_,n1), eye(k_), zeros(k_,n2)]] ;
GAM1   = [GAM1t; [zeros(k_,n1+n2),eye(k_)]] ;
PSI    = [PSIt; zeros(k_,l_)] ;
PPI    = [zeros(n1+n2,k_); eye(k_)] ;
C      = [Ct; zeros(k_,1)] ;

[S0, S1, S2, ~, ~, ~, ~, unq] = smats(C, GAM0, GAM1, PSI, PPI);
%keyboard
out.unq = unq ; 

if unq<1
    out.Q = [] ; 
    out.G = [] ;
    return ; 
end

out.Q = [S1(smats2bp_ord,smats2bp_ord), ...
    -S0(smats2bp_ord,1) ; zeros(1,n_-1), 1] ;
out.G = [S2(smats2bp_ord,:) ; zeros(1,l_)] ;
out.unq = unq ;

out.bp2smats_ord = bp2smats_ord ;
out.smats2bp_ord = smats2bp_ord ;
