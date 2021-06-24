function [out, Ct, Pt, Dt, Qt, Gt, U] = loglike(xx, T_d, T_f, SET, Zobs)
%function out = loglike(params, T_d, T_f, SET, data)
%
% Calculates the log-likelihood for Gaussian sequence given a mean and
% variance sequence.
%
% State space representation:
% 
% S_t = C(t) + P(t) S_{t-1} + D(t) n_t : State-state equation
% X_t = H S_t + v_t : observation equation
%
% Zobs is nobs x T_sample

n_       = SET.variable.n_ ; 
nobs     = SET.nobs ;
T_sample = SET.ss; 
SC.t_d   = T_d ;
SC.t_f   = T_f ;

%% Structure 1: Re-solve with xx parameters

SET.M_.params(SET.EST.params_to_estimate_idx) = xx ;

% when estimating mu, or inflation targets, these lines are needed
%mu    = xx(find(contains(SET.EST.params_to_estimate,'mu'),1)) ;
%pit_d = xx(find(contains(SET.EST.params_to_estimate,'pit_d'),1)) ;
%pit_f = xx(find(contains(SET.EST.params_to_estimate,'pit_f'),1)) ;

mu      = SET.M_.params(SET.param_names.mu) ; 
pit_d   = SET.M_.params(SET.param_names.pit_d) ; 
pit_f   = SET.M_.params(SET.param_names.pit_f) ; 

bet = SET.M_.params(SET.param_names.bet) ;

SET.param_zlb_f     = log(pit_f) + log(mu) - log(bet) - 0.00125/4 ;
SET.param_zlb_d     = log(pit_d) + log(mu) - log(bet) - 0.0025/4;
SET.param_zlb_d_1pc = log(pit_d) + log(mu) - log(bet) - 0.01/4 ;

dyn_in.M_       = SET.M_ ;
dyn_in.oo_      = SET.oo_ ;
dyn_in.options_ = SET.options_ ;
dyn_in.solve    = 1 ;

mats = resolve(SET,dyn_in) ; % Construct structural matrices: A0, A1, etc

if isempty(mats.mats.Q) ; out = -10e12 ; return ; end

SET.mats.Q = mats.mats.Q ;
SET.mats.G = mats.mats.G ;

P = mats.mats.Q(1:end-1,1:end-1,:) ;
D = mats.mats.G(1:end-1,:,:) ;    

Omega  = eye(SET.variable.l_) ; SET.Omega = Omega ;
DOmega = D(:,:,1)*Omega*D(:,:,1)' ;
SET.Sigma = dlyap_doubling(P,DOmega) ;

%% Estimated Durations

SET.mat_init     = mats.mat_init ;
SET.mat_i_f_zlb  = mats.mat_i_f_zlb ;
SET.mat_i_d_zlb  = mats.mat_i_d_zlb ;
SET.mat_i_df_zlb = mats.mat_i_df_zlb ;
SET.mat_fin      = mats.mat_fin ;

mats = tv_mats(SET,SC) ; % Constructs Qt and Gt

%% Extract Matrices for Kalman filter

Qt = mats.Qt ;
Gt = mats.Gt ;

Ct = Qt(1:n_-1,end,:) ;
Pt = Qt(1:n_-1,1:n_-1,:) ;
Dt = Gt(1:n_-1,:,:) ;

H = zeros(nobs,n_-1) ; 

for i=1:nobs
    H(i,SET.EST.obs(i)) = 1 ;
end

%% Kalman Filter and Construct Likelihood

[U, St, flag] = kf_(SET, Ct, Pt, Dt, H, Zobs, SC) ;

if flag ; out=-1e12 ; return ; end

out = 0;

for t=1:T_sample
    
    hidx = SET.EST.hidx_no_r_f_no_r ;
    if     SC.t_f(t)==0 && SC.t_d(t)==0 ; hidx = SET.EST.hidx ;
    elseif SC.t_f(t)>0  && SC.t_d(t)==0 ; hidx = SET.EST.hidx_no_r_f ; 
    elseif SC.t_f(t)==0 && SC.t_d(t)>0  ; hidx = SET.EST.hidx_no_r ;
    end

    out = out - (nobs/2) * log(2*pi) - ...
        (1/2) * log(det(St(hidx,hidx,t))) - ...
        (1/2) * U(hidx,t)' / St(hidx,hidx,t) * U(hidx,t) ;

end

if abs(imag(out))>0
    out=-1e12;
    disp('imag ll') ;
end
