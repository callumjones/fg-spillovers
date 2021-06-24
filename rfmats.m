function out = rfmats(SET, z_d_i, z_d_f, z_f_i, z_f_f)
%function out = rfmats(SET, z_d_i, z_d_f, z_f_i, z_f_f)
%
%

n = SET.variable.n_ ;
l = SET.variable.l_ ;

Qf = SET.mats.Q ; 

%% Setup

if ~z_d_i && ~z_d_f
    z_d_i = z_f_f ;
    z_d_f = z_f_f-1 ;
end

Tbstar = max(z_d_f, z_f_f) ;

Qhat = zeros(n,n,Tbstar);
Ghat = zeros(n,l,Tbstar);

A0i = zeros(n,n,Tbstar) ;
A1i = zeros(n,n,Tbstar) ;
B0i = zeros(n,n,Tbstar) ;
D0i = zeros(n,l,Tbstar) ;
D2i = zeros(n,l,Tbstar) ;
Gammai = zeros(l,l,Tbstar) ;
Ai     = zeros(n,n,Tbstar) ;
A0invi = zeros(n,n,Tbstar) ;
Bi     = zeros(n,n,Tbstar) ;

%%

t_vec = 1:1:Tbstar ;
d_zlb = t_vec ;
f_zlb = t_vec ;

d_zlb([1:(z_d_i-1), (z_d_f+1):Tbstar])=0 ;
f_zlb([1:(z_f_i-1), (z_f_f+1):Tbstar])=0 ;

for t_ = 1:Tbstar

    if t_==d_zlb(t_) && t_~=f_zlb(t_)

        A0i(:,:,t_) = SET.mat_i_d_zlb.A0 ;
        A1i(:,:,t_) = SET.mat_i_d_zlb.A1 ;
        B0i(:,:,t_) = SET.mat_i_d_zlb.B0 ;
        D0i(:,:,t_) = SET.mat_i_d_zlb.D0 ;
        D2i(:,:,t_) = SET.mat_i_d_zlb.D2 ;
        Gammai(:,:,t_) = SET.mat_i_d_zlb.Gamma ;

        Ai(:,:,t_) = SET.mats.Ai_d_zlb ; %tmp_d * A1i(:,:,t_) ; 
        A0invi(:,:,t_) = SET.mats.A0invi_d_zlb ;%eye(n,n) / A0i(:,:,t_) ;
        Bi(:,:,t_) = SET.mats.Bi_d_zlb ; %tmp_d * B0i(:,:,t_) ;         
        
    elseif t_~=d_zlb(t_) && t_==f_zlb(t_)

        A0i(:,:,t_) = SET.mat_i_f_zlb.A0 ;
        A1i(:,:,t_) = SET.mat_i_f_zlb.A1 ;
        B0i(:,:,t_) = SET.mat_i_f_zlb.B0 ;
        D0i(:,:,t_) = SET.mat_i_f_zlb.D0 ;
        D2i(:,:,t_) = SET.mat_i_f_zlb.D2 ;
        Gammai(:,:,t_) = SET.mat_i_f_zlb.Gamma ;
        
        Ai(:,:,t_) = SET.mats.Ai_f_zlb ; %tmp_d * A1i(:,:,t_) ; 
        A0invi(:,:,t_) = SET.mats.A0invi_f_zlb ;%eye(n,n) / A0i(:,:,t_) ;
        Bi(:,:,t_) = SET.mats.Bi_f_zlb ; %tmp_d * B0i(:,:,t_) ;         
        
    elseif t_==d_zlb(t_) && t_==f_zlb(t_)

        A0i(:,:,t_) = SET.mat_i_df_zlb.A0 ;
        A1i(:,:,t_) = SET.mat_i_df_zlb.A1 ;
        B0i(:,:,t_) = SET.mat_i_df_zlb.B0 ;
        D0i(:,:,t_) = SET.mat_i_df_zlb.D0 ;
        D2i(:,:,t_) = SET.mat_i_df_zlb.D2 ;
        Gammai(:,:,t_) = SET.mat_i_df_zlb.Gamma ;

        Ai(:,:,t_) = SET.mats.Ai_df_zlb ; %tmp_d * A1i(:,:,t_) ; 
        A0invi(:,:,t_) = SET.mats.A0invi_df_zlb ;%eye(n,n) / A0i(:,:,t_) ;
        Bi(:,:,t_) = SET.mats.Bi_df_zlb ; %tmp_d * B0i(:,:,t_) ;         
        
    else 
        
        SET.mat_i = SET.mat_init ; 

        A0i(:,:,t_) = SET.mat_i.A0 ;
        A1i(:,:,t_) = SET.mat_i.A1 ;
        B0i(:,:,t_) = SET.mat_i.B0 ;
        D0i(:,:,t_) = SET.mat_i.D0 ;
        D2i(:,:,t_) = SET.mat_i.D2 ;
        Gammai(:,:,t_) = SET.mat_i.Gamma ;
        
        Ai(:,:,t_) = SET.mats.Ai_no_zlb ; %tmp_d * A1i(:,:,t_) ; 
        A0invi(:,:,t_) = SET.mats.A0invi_no_zlb ;%eye(n,n) / A0i(:,:,t_) ;
        Bi(:,:,t_) = SET.mats.Bi_no_zlb ; %tmp_d * B0i(:,:,t_) ;         
        
    end 

end

Qtil = (eye(n)-Bi(:,:,end)*Qf)\Ai(:,:,end) ;

for t = Tbstar-1:-1:1
    tmp_a       = eye(n)-Bi(:,:,t)*Qtil ;
    tmp_b       = eye(n,n) / (A0i(:,:,t)-B0i(:,:,t)*Qtil) ;
    Qtil        = tmp_a \ Ai(:,:,t) ;
    Qhat(:,:,t) = tmp_b * A1i(:,:,t);
    Ghat(:,:,t) = tmp_b * D0i(:,:,t);
end

tmp_c = eye(n,n) / (A0i(:,:,end)-B0i(:,:,end)*Qf) ;

Qhat(:,:,Tbstar) = tmp_c * A1i(:,:,end);
Ghat(:,:,Tbstar) = tmp_c * D0i(:,:,end);

%%

out.Qhat = Qhat ;
out.Ghat = Ghat ;