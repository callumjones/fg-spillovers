function [x_T, n_T] = ks_(SET, Ct, Pt, Dt, H, data, SC)
%function [mean, var] = ks(SET, Ct, Pt, Dt, H, V, data, SC) 
%
% Kalman smoother
%

n_   = SET.variable.n_ ;
ss   = SET.ss ;
nobs = SET.nobs ;

%% Initialize

Omega  = SET.Omega ;
DOmega = Dt(:,:,1)*Omega*Dt(:,:,1)' ;

Sigma = dlyap_doubling(Pt(:,:,1),DOmega) ;

S     = zeros(nobs,nobs,ss) ;
y_til = zeros(nobs,ss) ; 
x     = zeros(n_-1,ss) ;
P     = zeros(n_-1,n_-1,ss) ;
K     = zeros(n_-1,nobs,ss) ;

Ph = Sigma ;
xh = ((eye(n_-1) - Pt(:,:,1))) \ Ct(:,:,1)  ;

P(:,:,1) = Ph ;
x(:,1)   = xh ;

%% Obtain Kalman filter matrices with time T information

for t=1:ss
    hidx = SET.EST.hidx_no_r_f_no_r ;
    if     SC.t_f(t)==0 && SC.t_d(t)==0 ; hidx = SET.EST.hidx ;
    elseif SC.t_f(t)>0  && SC.t_d(t)==0 ; hidx = SET.EST.hidx_no_r_f ; 
    elseif SC.t_f(t)==0 && SC.t_d(t)>0  ; hidx = SET.EST.hidx_no_r ;
    end

    % Update
    y_til(hidx,t) = data(hidx,t) - H(hidx,:)*x(:,t) ;
    S(hidx,hidx,t) = H(hidx,:)*P(:,:,t)*H(hidx,:)' ;
    K(:,hidx,t) = P(:,:,t)*H(hidx,:)' / S(hidx,hidx,t) ;
    xh = x(:,t) + K(:,hidx,t)*y_til(hidx,t);
    Ph = (eye(n_-1) - K(:,hidx,t)*H(hidx,:))*P(:,:,t) ; 
    
	% Predict
    if t<ss
        x(:,t+1) = Ct(:,:,t+1)+Pt(:,:,t+1)*xh ; 
        P(:,:,t+1) = Pt(:,:,t+1)*Ph*Pt(:,:,t+1)' + ...
            Dt(:,:,t+1)*Omega*Dt(:,:,t+1)' ;
    end
end

%% Kalman smoother with time T information

x_T(:,ss) = x(:,ss) ;
r_T(:,ss+1) = zeros(n_-1,1) ;

for t=ss:-1:1
    hidx = SET.EST.hidx_no_r_f_no_r ;
    if     SC.t_f(t)==0 && SC.t_d(t)==0 ; hidx = SET.EST.hidx ;
    elseif SC.t_f(t)>0  && SC.t_d(t)==0 ; hidx = SET.EST.hidx_no_r_f ; 
    elseif SC.t_f(t)==0 && SC.t_d(t)>0  ; hidx = SET.EST.hidx_no_r ;
    end

    r_T(:,t) = H(hidx,:)' / S(hidx,hidx,t) ...
        * (data(hidx,t) - H(hidx,:)*x(:,t) ) + ...
        (eye(n_-1) - K(:,hidx,t)*H(hidx,:))' * Pt(:,:,t)' * r_T(:,t+1) ;

    x_T(:,t) = x(:,t) + P(:,:,t) * r_T(:,t) ;
    n_T(:,t) = Dt(:,:,t)' * r_T(:,t) ;
end
