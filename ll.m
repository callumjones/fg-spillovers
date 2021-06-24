function [out, U] = ll(xx, T_d, T_f, SET, data)
%function out = ll(params, SET, data)
%
% Calculates the log-likelihood for Gaussian sequence given a mean and
% variance sequence.
%
% State space representation:
% 
% S_t = C(t) + P(t) S_{t-1} + D(t) n_t : State-state equation
% X_t = H S_t + v_t : observation equation

n    = SET.variable.n ; 
nobs = SET.nobs ;

%%

%xx(21)=10;

SET.param = xx ; 
SC.t_d    = T_d ;
SC.t_f    = T_f ;

%SET.param_calib.mu  = xx(30) ;
%SET.param_pit_d     = xx(31) ;
%SET.param_pit_f     = xx(32) ;

SET.param_zlb_f         = log(SET.param_pit_f) + ...
    log(SET.param_calib.mu) - log(SET.param_calib.beta) - 0.00125/4 ;
SET.param_zlb_d         = log(SET.param_pit_d) + ...
    log(SET.param_calib.mu) - log(SET.param_calib.beta) - 0.0025/4;
SET.param_zlb_d_1pc     = log(SET.param_pit_d) + ...
    log(SET.param_calib.mu) - log(SET.param_calib.beta) - 0.01/4 ;

switch SET.kf
    case 0
        mats = tv_mats(SET,SC) ;
    case 1
        FC.t_d = [0 SC.t_d] ;
        FC.t_f = [0 SC.t_f] ;

        FC.t_d(FC.t_d>0) = FC.t_d(FC.t_d>0)-1 ;
        FC.t_f(FC.t_f>0) = FC.t_f(FC.t_f>0)-1 ;

        mats = tv_mats(SET,FC) ;
end

%keyboard

if mats.unq==0 ; out=-1e8 ; return ; end

Qt = mats.Qt ;
Gt = mats.Gt ;

Ct = Qt(1:n-1,end,:);
Pt = Qt(1:n-1,1:n-1,:);
Dt = Gt(1:n-1,:,:);

H = zeros(nobs,n-1) ; 

for i=1:nobs
    H(i,SET.EST.obs_(i)) = 1 ;
end

[U, St, flag] = kf_(SET, Ct, Pt, Dt, H, data, SC) ;

if flag ; out=-1e12 ; return ; end

out = 0;

if SET.init_Sigma==0 
    st=2 ;
elseif SET.init_Sigma==1
    st=SET.EST.st ;
end

for t=st:length(data)
    if SC.t_f(t)==0&&SC.t_d(t)==0 ; hidx = SET.EST.hidx ;
    elseif SC.t_f(t)>0&&SC.t_d(t)==0; hidx = SET.EST.hidx_no_r_f ; 
    elseif SC.t_f(t)==0&&SC.t_d(t)>0; hidx = SET.EST.hidx_no_r ;
    else hidx = SET.EST.hidx_no_r_f_no_r ;
    end
    
    out = out - (nobs/2) * log(2*pi) - ...
        (1/2) * log(det(St(hidx,hidx,t))) - ...
        (1/2) * U(hidx,t)' / St(hidx,hidx,t) * U(hidx,t) ;
end

if abs(imag(out))>0
    out=-1e12;
    disp('imag ll') ;
end
