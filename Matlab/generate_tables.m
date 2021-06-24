% 
% Jones, Kulish, Rees
% International Spillovers of Forward Guidance Shocks
% Journal of Applied Econometrics, 2021
%

clear all
clc

%%

path(pathdef) ;

warning('Off','all') ;

addpath('./routines')
addpath('./routines/solmat')
addpath(genpath('./routines/VAR-Toolbox/v3dot0/'))
addpath('./output') 
addpath('./input')

%% Load posterior chain and find mode

load mhall_05-Mar-2020_estimatechis.mat

params_x_   = [] ;
params_T_d_ = [] ;
params_T_f_ = [] ;

burn = 50000 ;

for runs=1:maxproc
    params_x_ = [params_x_ ; params_x(burn:end,:,runs)] ; 
    params_T_d_ = [params_T_d_ ; params_T_d(burn:end,:,runs)] ;
    params_T_f_ = [params_T_f_ ; params_T_f(burn:end,:,runs)] ; 
end

for ii=1:length(params_x_(1,:))
    %[N,X] = hist(params_x_(:,ii),5000) ;
    [a,b]=ksdensity(params_x_(:,ii),'NumPoints',200);
    %tmp = X(N==max(N)) ;
    x_mode_chain(ii) = b(a==max(a)) ;
end

Tbs_out_d_med = round(mean(params_T_d_)) ;
Tbs_out_f_med = round(mean(params_T_f_)) ;

pick_draw = 125000; 
x_draw = params_x_(pick_draw,:); 
T_d_draw = params_T_d_(pick_draw,:); 
T_f_draw = params_T_f_(pick_draw,:); 

%save ./output/x_T_d_T_f x_mode_chain Tbs_out_d_med Tbs_out_f_med

%save ./output/onedraw x_draw T_d_draw 

%% Print Table

table = [] ;

%priorx = draw_from_prior(SET, params_x(end,:,1)) ; save ./output/prior_sav priorx;
load prior_sav % if done before

table = [table prctile(priorx,[50 10 90])] ;
table = [table' x_mode_chain' prctile(params_x_,[50 10 90])'] ;

disp('Prior Median, 10pc, 90pc, and Posterior Mode, Median, 10pc, 90pc, of parameters') ;
latex(table, '%.2f', 'nomath' ) ;

%% Print Variance Decomposition Table

clear all

setup_model

%     'eps_z     '  3
%     'eps_xi_f  '  1
%     'eps_r_f   '  2
%     'eps_g_f   '  6
%     'eps_xi_p_f'  5
%     'eps_xi_w_f'  7
%     'eps_xi    '  8
%     'eps_r     '  9
%     'eps_g     '  11
%     'eps_rp    '  10
%     'eps_xi_H  '  15
%     'eps_xi_w  '  16
%     'eps_xi_X  '  13 
%     'eps_xi_F  '  14

table = [] ;
table = oo_.variance_decomposition(:,[3 1 2 6 5 7 8 9 11 10 15 16 13 14]) ;
latex(table, '%.1f', 'nomath' ) ;

% Ordering of variables
% r_f_1_obs r_f_8_obs dy_f_obs dc_f_obs infl_f_obs winf_f_obs r_1_obs r_8_obs dy_obs dc_obs infl_obs winf_obs dy_F_obs dx_obs de_obs
