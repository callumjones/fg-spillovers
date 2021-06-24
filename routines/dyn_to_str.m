function out = dyn_to_str(dyn_in)
%{

The linearized model is:
 A x(t) = B x(t-1) + D E_t x(t+1) + E e(t) 

The reduced form model is:
 x(t) = F x(t-1) + G e(t) 

The constant is included in the structural matrices.

The approach will be:
 1. Run the Dynare model
 2. Extract the Jacobian from the 'dynamic model' specification
 3. Construct the structural matrices from the Jacobian
 4. Use those structural matrices with the ZLB algorithm

Callum Jones
callum.jones@nyu.edu 

%}
%%
%{
Why the Jacobian? From the Dynare manual: 

When computing the Jacobian of the dynamic model, the order of the
endogenous variables in the columns is stored in M_.lead_lag_incidence. The
rows of this matrix represent time periods: the first row denotes a lagged
(time t-1) variable, the second row a contemporaneous (time t) variable,
and the third row a leaded (time t+1) variable. The columns of the matrix
represent the endogenous variables in their order of declaration. A zero in
the matrix means that this endogenous does not appear in the model in this
time period. The value in the M_.lead_lag_incidence matrix corresponds to
the column of that variable in the Jacobian of the dynamic model. Example:
Let the second declared variable be c and the (3,2) entry of
M_.lead_lag_incidence be 15. Then the 15th column of the Jacobian is the
derivative with respect to c(+1).

Comparison to Dynare output: 

F(:,oo_.dr.state_var) is equivalent to oo_.dr.ghx(oo_.dr.inv_order_var,:)
G is equivalent to G2 where
- G2=oo_.dr.ghu*sqrt(M_.Sigma_e) ;
- G2(oo_.dr.inv_order_var,:) ;

%}
%%

M_         = dyn_in.M_ ;
oo_        = dyn_in.oo_ ;
options_   = dyn_in.options_ ;
model_name = cellstr(M_.fname) ;

%%

params = M_.params ; %#ok<NASGU>
M_ind  = M_.lead_lag_incidence ;

n_ = M_.endo_nbr ;       % Number of endogenous variables
l_ = M_.exo_nbr ;        % Number of exogenous variables
n2 = sum(M_ind(3,:)>0) ; % Number of expectational variables
n1 = n_-n2 ;             % Number of non-expectational variables

%% Compute Jacobian matrix

% tmp does some indexing needed to work out the state vector at
% which the Jacobian is evaluated.
tmp = repmat(1:n_,[3,1]) ;
tmp = tmp.*(M_.lead_lag_incidence>0) ;
tmp = tmp' ;
tmp = tmp(:) ;
tmp = tmp(tmp>0) ;

% if linearizing around another point, replace ss with that point
if dyn_in.linearize_around_diff_y
    ss = dyn_in.dyn_linearize_point ;
else
   	[ss, params, info] = steady_(M_,options_,oo_);
end

%if abs(ss(12)-log(M_.params(8)))>1e-8 ; disp('indeterminate') ; return ; end
%keyboard

y = ss(tmp) ; %#ok<NASGU>
x = zeros(1,l_) ; %#ok<NASGU>

eval(['[~, jac] = ', model_name{1}, '.dynamic(y, x, params, [], 1) ;']) ;

%% Extract structural matrices from Jacobian

FF = zeros(n_,n_) ;
GG = zeros(n_,n_) ;
HH = zeros(n_,n_) ;
MM = zeros(n_,l_) ;  % #ok<NASGU>

yind = 1:n_ ;

tmp1  = M_ind(1,:) ; % The 1st row will define the t-1 matrix
yind1 = yind(tmp1>0) ;
tmp1(tmp1==0)=[] ;

tmp2  = M_ind(2,:) ; % The 2nd row will define the t matrix
yind2 = yind(tmp2>0) ;
tmp2(tmp2==0)=[] ;

tmp3  = M_ind(3,:) ; % The 3rd row will define the t+1 matrix
yind3 = yind(tmp3>0) ;
tmp3(tmp3==0)=[] ;

GG(:,yind2) = jac(:,tmp2) ; %#ok<NODEF>
HH(:,yind1) = -jac(:,tmp1) ;
FF(:,yind3) = -jac(:,tmp3) ;
MM          = -jac(:,max(M_ind(:))+1:end)*sqrt(M_.Sigma_e) ;

%% Add constant

n_ = n_+1 ;

C = (GG-HH-FF)*ss ;
GG(1:n_-1,n_) = -C(1:n_-1) ;

GG(n_,n_) = 1 ; 
HH(n_,n_) = 1 ; 
FF(n_,n_) = 0 ;
MM(n_,:)  = 0 ;

%% Solve the system

A    = GG ;
B    = HH ;
D    = FF ;
E    = MM ;
D2    = zeros(n_,l_) ;
Gamma = zeros(l_,l_) ;

mats = bp2smats(M_,A,B,D,E) ;

F = mats.Q ;
G = mats.G ;

%[F, G, ~, ~, ~] = ...
%    smatsbp(A, B, D, E, D2, Gamma) ;

%% Outputs

out.mats.A = A ;
out.mats.B = B ;
out.mats.D = D ;
out.mats.E = E ;
out.mats.D2 = D2 ;
out.mats.Gamma = Gamma ;
out.mats.Q = F ;
out.mats.G = G ;
out.mats.unq = mats.unq ;
