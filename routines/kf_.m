function [y_til, S, flag, Sigma] = kf_(SET, Ct, Pt, Dt, H, data, SC)
%function [mean, var] = kf(SET, Ct, Pt, Dt, H, data, SC) 
% 
% Kalman Filter
%

n_   = SET.variable.n_ ;
ss   = SET.ss ;
nobs = SET.nobs ;
flag = 0 ;

%% Initialize

Sigma = SET.Sigma ; 
Omega = SET.Omega ;

S     = zeros(nobs,nobs,ss) ;
y_til = zeros(nobs,ss) ; 
x     = zeros(n_-1,nobs) ;
P     = zeros(n_-1,n_-1,nobs) ;
K     = zeros(n_-1,nobs,ss) ;

Ph = Sigma ;
xh = ((eye(n_-1) - Pt(:,:,1))) \ Ct(:,:,1)  ;

P(:,:,1) = Ph ;
x(:,1)   = xh ;

%% Kalman Filter

for t=1:ss
    hidx = SET.EST.hidx_no_r_f_no_r ;
    if     SC.t_f(t)==0  && SC.t_d(t)==0 ; hidx = SET.EST.hidx ;
     elseif SC.t_f(t)>0  && SC.t_d(t)==0 ; hidx = SET.EST.hidx_no_r_f ; 
     elseif SC.t_f(t)==0 && SC.t_d(t)>0  ; hidx = SET.EST.hidx_no_r ;
    end

    % Update
    y_til(hidx,t) = data(hidx,t) - H(hidx,:)*x(:,t) ;
    
    S(hidx,hidx,t) = H(hidx,:)*P(:,:,t)*H(hidx,:)' ; 
    K(:,hidx,t) = P(:,:,t)*H(hidx,:)' / S(hidx,hidx,t)  ;
    xh = x(:,t) + K(:,hidx,t)*y_til(hidx,t);
    Ph = (eye(n_-1) - K(:,hidx,t)*H(hidx,:))*P(:,:,t) ; 
    
	% Predict
    if t<ss 
        x(:,t+1)   = Ct(:,:,t+1)+Pt(:,:,t+1)*xh ;
        P(:,:,t+1) = Pt(:,:,t+1)*Ph*Pt(:,:,t+1)' + ...
                     Dt(:,:,t+1)*Omega*Dt(:,:,t+1)' ;
    end
end
