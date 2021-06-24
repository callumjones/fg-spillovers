function out = resolve(SET,dyn_in)
%%

set_params

dyn_in.M_.params = out ;

out = [] ;

%dyn_in.dyn_linearize_point = xsteady' ;
dyn_in.linearize_around_diff_y = 0 ;

%%

mats_out = dyn_to_str(dyn_in) ;

if dyn_in.solve
out.mats.Q = mats_out.mats.Q ;
out.mats.G = mats_out.mats.G ; 
out.zlb.Qf = out.mats.Q ; 
end

out.str_mats.A0 = mats_out.mats.A ;
out.str_mats.A1 = mats_out.mats.B ;
out.str_mats.B0 = mats_out.mats.D ;
out.str_mats.D0 = mats_out.mats.E ;
out.str_mats.D2 = mats_out.mats.D2 ;
out.str_mats.Gamma = mats_out.mats.Gamma ;

% ZLB only in large country
out.mat_i_f_zlb = out.str_mats ;

out.mat_i_f_zlb.A0(SET.f_tr,:) = 0 ;
out.mat_i_f_zlb.A1(SET.f_tr,:) = 0 ;
out.mat_i_f_zlb.B0(SET.f_tr,:) = 0 ;
out.mat_i_f_zlb.D0(SET.f_tr,:) = 0 ;
out.mat_i_f_zlb.A0(SET.f_tr,SET.variable.r_f_1) = 1 ;
out.mat_i_f_zlb.A0(SET.f_tr,end) = SET.param_zlb_f ;

% ZLB only in small country
out.mat_i_d_zlb = out.str_mats ;

out.mat_i_d_zlb.A0(SET.d_tr,:) = 0 ;
out.mat_i_d_zlb.A1(SET.d_tr,:) = 0 ;
out.mat_i_d_zlb.D0(SET.d_tr,:) = 0 ;
out.mat_i_d_zlb.B0(SET.d_tr,:) = 0 ;
out.mat_i_d_zlb.A0(SET.d_tr,SET.variable.r_1) = 1 ;
out.mat_i_d_zlb.A0(SET.d_tr,end) = SET.param_zlb_d ;

% ZLB in both large and small countries
out.mat_i_df_zlb = out.str_mats ;

out.mat_i_df_zlb.A0(SET.d_tr,:) = 0 ;
out.mat_i_df_zlb.A1(SET.d_tr,:) = 0 ;
out.mat_i_df_zlb.D0(SET.d_tr,:) = 0 ;
out.mat_i_df_zlb.B0(SET.d_tr,:) = 0 ;
out.mat_i_df_zlb.A0(SET.d_tr,SET.variable.r_1) = 1 ;
out.mat_i_df_zlb.A0(SET.d_tr,end) = SET.param_zlb_d ;

out.mat_i_df_zlb.A0(SET.f_tr,:) = 0 ;
out.mat_i_df_zlb.A1(SET.f_tr,:) = 0 ;
out.mat_i_df_zlb.D0(SET.f_tr,:) = 0 ;
out.mat_i_df_zlb.B0(SET.f_tr,:) = 0 ;
out.mat_i_df_zlb.A0(SET.f_tr,SET.variable.r_f_1) = 1 ;
out.mat_i_df_zlb.A0(SET.f_tr,end) = SET.param_zlb_f ;

out.mat_init    = out.str_mats ; % Structural matrices not at ZLB
out.mat_fin     = out.str_mats ; % Structural matrices after ZLB  
