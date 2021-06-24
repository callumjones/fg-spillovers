% Simulates the model under the ZLB
%
%
%%

clear all
clc

%%

path(pathdef) ;

warning('Off','all') ;

addpath('./routines')
addpath('./routines/solmat')

setup_model

%%

SET.EST.params_to_estimate = ...
    {'h_f', ...
    'theta_p_f', 'theta_w_f', ...
     'rho_r_f', 'phi_pi_f', 'phi_g_f', 'phi_y_f', ...
     'rho_xi_f', 'rho_g_f', 'rho_xi_p_f', ...
     'rho_xi_w_f','rho_tp_f_8',...
     'sig_z','sig_r_f', 'sig_xi_f', 'sig_g_f' , 'sig_xi_p_f', ...
     'sig_xi_w_f', 'sig_eps_r_f_8', 'r_f_8_const', ...
     'h','tau',...
     'theta_p','theta_w', ...
     'theta_x', 'theta_F', ...
        'rho_r', 'phi_pi', 'phi_g', 'phi_y', ...
        'rho_rp','rho_xi','rho_g',...
        'rho_xi_H','rho_xi_w','rho_xi_X','rho_xi_F','rho_tp_8', ...
        'sig_r','sig_rp','sig_g',...
        'sig_xi','sig_xi_H','sig_xi_w', 'sig_xi_X','sig_xi_F',...
        'sig_eps_r_8', 'r_8_const'  ...
        'chi_p_f', 'chi_w_f', 'chi_p', 'chi_w'} ;


SET.EST.params_to_estimate_idx = [] ; 
 
for qq=1:length(SET.EST.params_to_estimate)
    tmp = find(strcmp(cellstr(dyn_in.M_.param_names),SET.EST.params_to_estimate{qq})) ;
    SET.EST.params_to_estimate_idx(qq) = tmp ;
end

SET.EST.varobs = { ...
    'r_f_1_obs', 'r_f_8_obs', 'dy_f_obs','dc_f_obs', 'infl_f_obs','winf_f_obs', ...
    'r_1_obs', 'r_8_obs', 'dy_obs','dc_obs','infl_obs', 'dy_F_obs', 'dx_obs', 'de_obs', 'winf_obs'} ;

SET.EST.obs = [] ; 
 
for qq=1:length(SET.EST.varobs)
    tmp = find(strcmp(cellstr(SET.M_.endo_names),SET.EST.varobs{qq})) ;
    SET.EST.obs(qq) = tmp ;
end

dyn_in.M_  = M_ ;
dyn_in.oo_ = oo_ ;
dyn_in.solve = 1 ;
dyn_in.linearize_around_diff_y = 0 ;

load x_T_d_T_f.mat

% sig_a_f is in 10
% sig_xi_f is in 11
% sig_xi_p_f is in 13

%x_mode_fmin([11])=0.01;
%load hessmat
xx_from_chain = x_mode_chain;

% xx_from_chain(17) =  xx_from_chain(17) ;
% xx_from_chain(42) =  xx_from_chain(42) ;
% xx_from_chain(43) =  xx_from_chain(43) ;

%dyn_in.M_.params(SET.EST.params_to_estimate_idx) = x_mode_fmin ;
dyn_in.M_.params(SET.EST.params_to_estimate_idx) = xx_from_chain ;


SET.M_ = dyn_in.M_ ;
SET.oo_ = dyn_in.oo_ ;
SET.options_ = options_ ;

out = dyn_to_str(dyn_in) ;

SET.mats.Q = out.mats.Q ; Q = SET.mats.Q ;
SET.mats.G = out.mats.G ; G = SET.mats.G ;


%% Constructs a stochastic simulation

k = 2 ; % State of random number generator

randn('state',k) ;

SET.horizon = 250 ; % Length of the simulation
Q = SET.mats.Q ;
G = SET.mats.G ;

et = 1*(eye(SET.variable.l_,SET.variable.l_).^.5) * ...
    randn(SET.variable.l_,SET.horizon) ; % Draw shocks

% Compute steady-state from structural matrices
SET.y_ss = (eye(SET.variable.n_-1) - Q(1:SET.variable.n_-1,1:SET.variable.n_-1)) ...
    \ Q(1:SET.variable.n_-1,end)  ;
SET.y_ss = [SET.y_ss ; 1] ; % Constant

SET.EST.zlb_at_1pc = zeros(1,SET.horizon) ;

%% Run simulation

z_f_T    = zeros(1,SET.horizon) ;
z_d_T    = zeros(1,SET.horizon) ;
zlb_f    = zeros(1,SET.horizon) ;
zlb_d    = zeros(1,SET.horizon) ;
y_no_zlb = SET.y_ss ;
y_zlb    = SET.y_ss ;

for t = 2:SET.horizon
    
    y_no_zlb(:,t) = Q*y_no_zlb(:,t-1) + G*et(:,t) ;
    
    [y_zlb(:,t), zlb_f(t), zlb_d(t), z_f_T(t), z_d_T(t)] = ... 
        zlb(SET, y_zlb(:,t-1), et(:,t)) ;
    
end

%% 

load_data

figure ;
for ii=1:4 ; % SET.nobs 
    subplot(2,2,ii) ; hold on ;
    plot(data_dates,Zobs(ii,:)) ;
    plot(data_dates(1:111),y_no_zlb(SET.EST.obs(ii),1:111));
    title(SET.EST.varobs{ii}) ;
    xlim([data_dates(1) data_dates(end)]) ;
end
print -depsc ./sim_v_data1.eps
%close


figure ;
for ii=5:8 ; % SET.nobs 
    subplot(2,2,ii-4) ; hold on ;
    plot(data_dates,Zobs(ii,:)) ;
    plot(data_dates(1:111),y_no_zlb(SET.EST.obs(ii),1:111));
    title(SET.EST.varobs{ii}) ;
    xlim([data_dates(1) data_dates(end)]) ;
end
print -depsc ./sim_v_data2.eps


figure ;
for ii=9:12 ; % SET.nobs 
    subplot(2,2,ii-8) ; hold on ;
    plot(data_dates,Zobs(ii,:)) ;
    plot(data_dates(1:111),y_no_zlb(SET.EST.obs(ii),1:111));
    title(SET.EST.varobs{ii}) ;
    xlim([data_dates(1) data_dates(end)]) ;
end
print -depsc ./sim_v_data3.eps

figure ;
for ii=13:15 ; % SET.nobs 
    subplot(2,2,ii-12) ; hold on ;
    plot(data_dates,Zobs(ii,:)) ;
    plot(data_dates(1:111),y_no_zlb(SET.EST.obs(ii),1:111));
    title(SET.EST.varobs{ii}) ;
    xlim([data_dates(1) data_dates(end)]) ;
end
print -depsc ./sim_v_data4.eps
