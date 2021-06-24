function out = tv_mats(SET,SC)
%function out = tv_mats(SET,SC)
%
%

ss = SET.ss ;
n  = SET.variable.n_ ;

inv_A0i         = eye(n,n) / SET.mat_init.A0 ;
inv_A0i_d_zlb   = eye(n,n) / SET.mat_i_d_zlb.A0 ;
inv_A0i_f_zlb   = eye(n,n) / SET.mat_i_f_zlb.A0 ;
inv_A0i_df_zlb  = eye(n,n) / SET.mat_i_df_zlb.A0 ;

SET.mats.Ai_no_zlb     = inv_A0i * SET.mat_init.A1 ;
SET.mats.A0invi_no_zlb = eye(n,n) / SET.mat_init.A0 ;
SET.mats.Bi_no_zlb     = inv_A0i * SET.mat_init.B0 ;

SET.mats.Ai_d_zlb     = inv_A0i_d_zlb * SET.mat_i_d_zlb.A1 ;
SET.mats.A0invi_d_zlb = eye(n,n) / SET.mat_i_d_zlb.A0 ;
SET.mats.Bi_d_zlb     = inv_A0i_d_zlb * SET.mat_i_d_zlb.B0 ;

SET.mats.Ai_f_zlb     = inv_A0i_f_zlb * SET.mat_i_f_zlb.A1 ;
SET.mats.A0invi_f_zlb = eye(n,n) / SET.mat_i_f_zlb.A0 ;
SET.mats.Bi_f_zlb     = inv_A0i_f_zlb * SET.mat_i_f_zlb.B0 ;

SET.mats.Ai_df_zlb     = inv_A0i_df_zlb * SET.mat_i_df_zlb.A1 ;
SET.mats.A0invi_df_zlb = eye(n,n) / SET.mat_i_df_zlb.A0 ;
SET.mats.Bi_df_zlb     = inv_A0i_df_zlb * SET.mat_i_df_zlb.B0 ;

for t_ = 1:ss
    
    % If FG at 1 pc in small country, change matrices
    if SET.EST.zlb_at_1pc(t_)>0 && SET.EST.zlb_at_1pc(t_-1)==0 
        SET.mat_i_d_zlb.A0(SET.d_tr,end) = SET.param_zlb_d_1pc ;
        
        inv_A0i_d_zlb         = eye(n,n) / SET.mat_i_d_zlb.A0 ;
        SET.mats.Ai_d_zlb     = inv_A0i_d_zlb * SET.mat_i_d_zlb.A1 ;
        SET.mats.A0invi_d_zlb = eye(n,n) / SET.mat_i_d_zlb.A0 ;
        SET.mats.Bi_d_zlb     = inv_A0i_d_zlb * SET.mat_i_d_zlb.B0 ;
        
        SET.mat_i_df_zlb.A0(SET.d_tr,end) = SET.param_zlb_d_1pc ;
        
        inv_A0i_df_zlb         = eye(n,n) / SET.mat_i_df_zlb.A0 ;
        SET.mats.Ai_df_zlb     = inv_A0i_df_zlb * SET.mat_i_df_zlb.A1 ;
        SET.mats.A0invi_df_zlb = eye(n,n) / SET.mat_i_df_zlb.A0 ;
        SET.mats.Bi_df_zlb     = inv_A0i_df_zlb * SET.mat_i_df_zlb.B0 ;  
    end
    
    if SC.t_f(t_)==0 && SC.t_d(t_)==0 % No ZLB
        
        out.Qt(:,:,t_) = SET.mats.Q ;
        out.Gt(:,:,t_) = SET.mats.G ;        
        
    elseif SC.t_f(t_)>0 && SC.t_d(t_)==0 % ZLB in large country
        
        z_f_i = 1 ;
        z_f_f = SC.t_f(t_) ;
        
        mats = rfmats(SET, 0, 0, z_f_i, z_f_f) ;
        
        out.Qt(:,:,t_) = mats.Qhat(:,:,1) ;
        out.Gt(:,:,t_) = mats.Ghat(:,:,1) ;
        
    elseif SC.t_f(t_)==0 && SC.t_d(t_)>0 % ZLB in small country
        
        z_d_i = 1 ;
        z_d_f = SC.t_d(t_) ;
        
        mats = rfmats(SET, z_d_i, z_d_f, 0, 0) ;
        
        out.Qt(:,:,t_) = mats.Qhat(:,:,1) ;
        out.Gt(:,:,t_) = mats.Ghat(:,:,1) ;
        
    else % ZLB in both large and small countries
        
        z_d_i = 1 ;
        z_d_f = SC.t_d(t_) ;
        z_f_i = 1 ;
        z_f_f = SC.t_f(t_) ; 
        
        mats = rfmats(SET, z_d_i, z_d_f, z_f_i, z_f_f) ;
        
        out.Qt(:,:,t_) = mats.Qhat(:,:,1) ;
        out.Gt(:,:,t_) = mats.Ghat(:,:,1) ;
        
    end
end
