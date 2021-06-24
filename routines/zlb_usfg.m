function [out, zlb_f, zlb_d, z_f_T, z_d_T] = zlb_usfg(SET, in, e_in, tt_in, FG)
%function [out, zlb_f, zlb_d, z_f_T, z_d_T] = zlb(SET, in, e_in, FG)
%
%

n = SET.variable.n_ ;

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

if SET.EST.zlb_at_1pc(tt_in)>0 %If FG at 1 pc
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

%%

z_f_i = FG.z_f_i ; 
z_f_f = FG.z_f_f ; 
z_d_i = FG.z_d_i ; 
z_d_f = FG.z_d_f ; 
            
if and(z_f_f==0,z_d_f==0) 
    Qtb = repmat(SET.mats.Q,[1 1 SET.maxchk+1]) ; 
    Gtb = repmat(SET.mats.G,[1 1 SET.maxchk+1]) ; 
else
    mats = rfmats(SET, z_d_i, z_d_f, z_f_i, z_f_f) ;
    
    Qtb = mats.Qhat ;
    Gtb = mats.Ghat ;
end

%% Setup

maxchk = SET.maxchk ;

ss = size(Qtb,3) ;

Qtb(:,:,(ss+1):(ss+SET.maxchk)) = repmat(SET.mats.Q,[1 1 SET.maxchk]) ;

i_f = SET.i_f ;
i_d = SET.i_d ;

%% Check for violation of ZLB

ychk = Qtb(:,:,1) * in + Gtb(:,:,1) * e_in ;

for t_ = 2:maxchk
    ychk(:,t_) = Qtb(:,:,t_) * ychk(:,t_-1) ;
end

if SET.EST.zlb_at_1pc(tt_in)>0 %If FG at 1 pc
    SET.param_zlb_d = SET.param_zlb_d_1pc ;
end  

chg_d = ychk(i_d,:)<-SET.param_zlb_d-1e-10 ;
chk_d = sum(chg_d) ;
chk   = chk_d ; %keyboard

zlb_f = 0 ; z_f_T = 0 ;

if ~chk
    switch SET.stoch
        case 1 ; out = ychk(:,1) ;
        case 0 ; out = ychk(:,1:maxchk) ;
    end
    zlb_f = 0 ; zlb_d = 0 ;
    z_f_T = 0 ; z_d_T = 0 ;
    return
end

%% Domestic block

z_d_i = 100 ;
z_d_f = 0 ;

max_t = max(z_d_f,z_f_f) ;

chg_d = ychk(i_d,:)<-SET.param_zlb_d-1e-10 ;
chk_d = sum(chg_d) ;

zlb_d = chk_d>0 ;

while chk_d
    
    for t_ = 1:length(chg_d)
        if chg_d(t_)==1
            
            z_d_i = min(z_d_i,t_) ;
            z_d_f = max(z_d_f,t_) ;
            max_t = max(max_t, z_d_f) ;
            
            mats = rfmats(SET, z_d_i, z_d_f, z_f_i, z_f_f) ;
            
            Qtb = mats.Qhat ;
            Gtb = mats.Ghat ;
             
            break
        end
    end
    
    ss = size(Qtb,3) ;

    Qtb(:,:,(ss+1):(ss+SET.maxchk)) = repmat(SET.mats.Q,[1 1 SET.maxchk]) ;

    ychk(:,1) = Qtb(:,:,1) * in + Gtb(:,:,1) * e_in ;
    
    for tt = 2:max(z_d_f,z_f_f)
        ychk(:,tt) = Qtb(:,:,tt) * ychk(:,tt-1) ;
    end
    
    for tt = (max(z_d_f,z_f_f)+1):maxchk
        ychk(:,tt) = Qtb(:,:,tt) * ychk(:,tt-1) ;
    end
    
chg_d = ychk(i_d,:)<-SET.param_zlb_d-1e-10 ;    
    
chk_d = sum(chg_d)>0 ;    
    
end

if zlb_d
    z_d_T = z_d_f - z_d_i + 1 ;
else
    z_d_T = 0 ;
end

%% Outputs

switch SET.stoch
    case 1 ; out = ychk(:,1) ; 
    case 0 ; out = ychk(:,1:maxchk) ;     
end