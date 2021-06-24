function out = rfmats(SET, z_f_i, z_f_f)
%function out = rfmats(SET, z_f_i, z_f_f)
%
%

n_ = SET.variable.n_ ;
l_ = SET.variable.l_ ;
Qf = SET.zlb.Qf ; 

%% Setup

Tbstar = z_f_f ;

if SET.tv_str_chg.flag
    cur_t  = SET.tv_str_chg.cur_t ;
    Tbstar = max(z_f_f,SET.tv_str_chg.Tbstar-cur_t+1) ;
    Qf     = SET.tv_str_chg.Qf ;
end

Qt = zeros(n_,n_,Tbstar);
Gt = zeros(n_,l_,Tbstar);

A0i = zeros(n_,n_,Tbstar) ;
A1i = zeros(n_,n_,Tbstar) ;
B0i = zeros(n_,n_,Tbstar) ;
D0i = zeros(n_,l_,Tbstar) ;

%%

t_vec = 1:1:Tbstar ;
f_zlb = t_vec ;

f_zlb([1:(z_f_i-1), (z_f_f+1):Tbstar])=0 ;

for t_ = 1:Tbstar
    mat_i   = SET.zlb.mat_init ; 
    mat_zlb = SET.zlb.mat_i_f_zlb ;
    
    if SET.tv_str_chg.flag 
        if length(SET.tv_str_chg.str_mats)<(t_+cur_t-1)
            mat_i   = SET.tv_str_chg.str_mats{end} ;
            mat_zlb = SET.tv_str_chg.mat_i_f_zlb{end} ;
        else
            mat_i   = SET.tv_str_chg.str_mats{t_+cur_t-1} ;
            mat_zlb = SET.tv_str_chg.mat_i_f_zlb{t_+cur_t-1} ;
        end
    end

    %save tmp
    A0i(:,:,t_) = mat_i.A0 ; 
    A1i(:,:,t_) = mat_i.A1 ;
    B0i(:,:,t_) = mat_i.B0 ;
    D0i(:,:,t_) = mat_i.D0 ;

    if t_==f_zlb(t_)
        A0i(:,:,t_) = mat_zlb.A0 ;
        A1i(:,:,t_) = mat_zlb.A1 ;
        B0i(:,:,t_) = mat_zlb.B0 ;
        D0i(:,:,t_) = mat_zlb.D0 ;
    end 
end

Qt(:,:,end) = (A0i(:,:,end)-B0i(:,:,end)*Qf)\A1i(:,:,end);
Gt(:,:,end) = (A0i(:,:,end)-B0i(:,:,end)*Qf)\D0i(:,:,end);

for tt = Tbstar-1:-1:1
    Qt(:,:,tt) = (A0i(:,:,tt)-B0i(:,:,tt)*Qt(:,:,tt+1))\A1i(:,:,tt);
    Gt(:,:,tt) = (A0i(:,:,tt)-B0i(:,:,tt)*Qt(:,:,tt+1))\D0i(:,:,tt);
end

%%

out.Qhat = Qt ;
out.Ghat = Gt ;