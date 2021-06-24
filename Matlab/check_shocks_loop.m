% 
% Jones, Kulish, Rees
% International Spillovers of Forward Guidance Shocks
% Journal of Applied Econometrics, 2021
%
% =========================================================================
% check_shocks_loop
%
% Loops over parameter draws. For each draw calculates the correlation
% between the empirical shocks.
% =========================================================================

%% 0. Preliminaries

clear ;
clc ;

rng('default') ;    % Reset random number generator for replicability

nloops = 1000 ;      % Number of parameter draws

%% 1. Set paths and load model and data and parameter draws

path(pathdef) ;

warning('Off','all') ;

addpath('./routines')
addpath('./routines/solmat')
addpath(genpath('./routines/VAR-Toolbox/v3dot0/'))
addpath('./output') 
addpath('./input')

setup_model

load_data

load mhall_05-Mar-2020_estimatechis.mat

%% 2. Set estimated parameters and durations

SET.EST.params_to_estimate = {...
    'h_f', ...
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

% Pick the parameters from posterior chain
% load params_from_chain
params_x_   = [] ;
params_T_d_ = [] ;
params_T_f_ = [] ;

burn = 50000 ;

for runs=1:maxproc
    params_x_ = [params_x_ ; params_x(burn:end,:,runs)] ; 
    params_T_d_ = [params_T_d_ ; params_T_d(burn:end,:,runs)] ;
    params_T_f_ = [params_T_f_ ; params_T_f(burn:end,:,runs)] ; 
end

%% 3. Select draws and calculate correlations

ndraws    = size(params_x_,1) ; 

idx_draws = ceil(ndraws .* rand(nloops,1)) ;

tp_f = strmatch('eps_r_f_8',M_.exo_names) ;
tp_d = strmatch('eps_r_8',M_.exo_names) ;

rp_f = strmatch('eps_xi_f',M_.exo_names,'exact') ;
rp_d = strmatch('eps_xi',M_.exo_names,'exact') ;

for j_ = 1 : nloops

pick_draw = idx_draws(j_) ; 


xx_from_chain = params_x_(pick_draw,:);
dyn_in.M_.params(SET.EST.params_to_estimate_idx) = xx_from_chain ;
dursest_f = params_T_f_(pick_draw,:) ;
dursest_d = params_T_d_(pick_draw,:) ;

mu      = dyn_in.M_.params(SET.param_names.mu) ; 
pit_d   = dyn_in.M_.params(SET.param_names.pit_d) ; 
pit_f   = dyn_in.M_.params(SET.param_names.pit_f) ; 
bet     = dyn_in.M_.params(SET.param_names.bet) ;

SET.param_zlb_f     = log(pit_f) + log(mu) - log(bet) - 0.00125/4 ;
SET.param_zlb_d     = log(pit_d) + log(mu) - log(bet) - 0.0025/4;
SET.param_zlb_d_1pc = log(pit_d) + log(mu) - log(bet) - 0.01/4 ;

SET.M_ = dyn_in.M_ ;
SET.oo_ = dyn_in.oo_ ;
SET.options_ = options_ ;

mats = resolve(SET,dyn_in) ; % Construct structural matrices: A0, A1, etc

SET.mat_init     = mats.mat_init ;
SET.mat_i_f_zlb  = mats.mat_i_f_zlb ;
SET.mat_i_d_zlb  = mats.mat_i_d_zlb ;
SET.mat_i_df_zlb = mats.mat_i_df_zlb ;
SET.mat_fin      = mats.mat_fin ;

SET.mats.Q = mats.mats.Q ; 
SET.mats.G = mats.mats.G ; 

%% Kalman smoother to extract shocks

SET.nobs = length(data(:,1)) ; 
SET.ss   = length(data(1,:)) ; 

SET.EST.zlb_at_1pc = zeros(1,SET.ss) ;
SET.EST.zlb_at_1pc(1,find(data_dates==2010.5):find(data_dates==2014.75))=1 ;

SET.EST.varobs = { ...
    'r_f_1_obs', 'r_f_8_obs', 'dy_f_obs','dc_f_obs', 'infl_f_obs','winf_f_obs', ...
    'r_1_obs', 'r_8_obs', 'dy_obs','dc_obs','infl_obs', 'dy_F_obs', 'dx_obs', 'de_obs', 'winf_obs'} ;

SET.EST.obs = [] ; 
 
for qq=1:length(SET.EST.varobs)
    tmp = find(strcmp(cellstr(SET.M_.endo_names),SET.EST.varobs{qq})) ;
    SET.EST.obs(qq) = tmp ;
end

ss      = SET.ss ;
n       = SET.variable.n_ ; 
nparam  = SET.nparam ;
nobs    = SET.nobs ;

SET.Omega = eye(SET.variable.l_) ;

SET.EST.hidx        = 1:length(SET.EST.obs) ;
SET.EST.hidx_no_r_f = SET.EST.hidx(2:end) ;
SET.EST.hidx_no_r   = SET.EST.hidx([1:6 8:end]) ;
SET.EST.hidx_no_r_f_no_r = SET.EST.hidx([2:6 8:end]) ;

zlb_t_d = zeros(1,SET.ss) ;
zlb_t_f = zeros(1,SET.ss) ;

% Canadian effective percent lower bound ends on 2015Q1
zlb_t_d(1,find(data_dates==2009.25):find(data_dates==2014.75)) = 1 ;
% US effective lower bound ends on 2015Q4
zlb_t_f(1,find(data_dates==2009):find(data_dates==2015.5))     = 1 ;
zlb_t = [zlb_t_d ; zlb_t_f] ;

SC.t_f = zeros(1,SET.ss) ;
SC.t_d = zeros(1,SET.ss) ;  

SC.t_f(zlb_t_f>0) = dursest_f ; 
SC.t_d(zlb_t_d>0) = dursest_d ; 

% Obtain matrices consistent with the durations and parameters
mats = tv_mats(SET,SC) ;

Qt = mats.Qt ;
Gt = mats.Gt ;

Ct = Qt(1:n-1,end,:);
Pt = Qt(1:n-1,1:n-1,:);
Dt = Gt(1:n-1,:,:);

H = zeros(nobs,n-1) ; 

for i=1:nobs
    H(i,SET.EST.obs(i)) = 1 ;
end
    
[x_T, n_T] = ks_(SET, Ct, Pt, Dt, H, data, SC) ;

x_T = [x_T ; ones(1,ss) ] ;
n_T_ = n_T ;

% Save correlations of shocks
ccloop(:,j_) = vec(corrcoef(n_T')) ; 

% Obtain endogenous durations and FG shocks
SET.y_ss_init = (eye(n-1) - Qt(1:n-1,1:n-1,1)) \ Qt(1:n-1,end,1)  ;
SET.y_ss_init = [SET.y_ss_init ; 1] ;

y_no_sto  = SET.y_ss_init ;
y_no_zlb  = x_T(:,1) ;
y_zlb     = x_T(:,1) ;
z_f_T     = zeros(1,SET.ss) ;
z_d_T     = zeros(1,SET.ss) ;

Qt_nozlb = repmat(Qt(:,:,1),[1 1 SET.ss+SET.maxchk]) ;
Gt_nozlb = repmat(Gt(:,:,1),[1 1 SET.ss+SET.maxchk]) ;

for t=2:SET.ss
    y_no_sto(:,t) = SET.mats.Q*y_no_sto(:,t-1) ;
    y_no_zlb(:,t) = SET.mats.Q*y_no_zlb(:,t-1) + SET.mats.G*n_T(:,t) ;
    
    [y_zlb(:,t), ~, ~, z_f_T(t), z_d_T(t)] = ...
        zlb(SET, y_zlb(:,t-1), n_T(:,t), Qt_nozlb, Gt_nozlb, t) ;
end

lb_T_vec(j_,:)     = z_f_T(zlb_t_f>0) ; % ZLB duration
lb_T_d_vec(j_,:)   = z_d_T(zlb_t_d>0) ; % ZLB duration

% Calculate FG shocks
fg_vec_f    = dursest_f - lb_T_vec(j_,:) ; 
fg_vec_d    = dursest_d - lb_T_d_vec(j_,:) ;

fg_shock_f = fg_vec_f(1) ; 
fg_shock_d = fg_vec_d(1) ; 

for t_ = 2 : size(fg_vec_f,2)  
    if fg_vec_f(t_-1) == 0
        fg_shock_f(t_,1) = fg_vec_f(t_) ;
    else
        fg_shock_f(t_,1) = fg_vec_f(t_) - fg_vec_f(t_-1) + 1 ;
    end
end

for t_ = 2 : size(fg_vec_d,2) 
    if fg_vec_d(t_-1) == 0
        fg_shock_d(t_,1) = fg_vec_d(t_) ;
    else
        fg_shock_d(t_,1) = fg_vec_d(t_) - fg_vec_d(t_-1) + 1 ;
    end

end

fg_shock_save_f(:,j_) = fg_shock_f;
fg_shock_save_d(:,j_) = fg_shock_d ;

% [fg_vec_f' fg_shock_save_f(:,j_)]
% [fg_vec_d' fg_shock_save_d(:,j_)]
% 
% pause

% Save term premium shocks for ZLB quarters
%durest_save_f(:,j_)  = dursest_f' ;
%durest_save_d(:,j_)  = dursest_d' ;
%durest_shock_save_f(:,j_)       = [dursest_f(1) ;[dursest_f(2:end)-dursest_f(1:end-1)]'] ;
%durest_shock_save_d(:,j_)       = [dursest_d(1) ;[dursest_d(2:end)-dursest_d(1:end-1)]'] ;

tpshock_save_f(:,j_) = n_T(tp_f,zlb_t_f>0)' ;
tpshock_save_d(:,j_) = n_T(tp_d,zlb_t_d>0)' ;

rpshock_save_f(:,j_) = n_T(rp_f,zlb_t_f>0)' ;
rpshock_save_d(:,j_) = n_T(rp_d,zlb_t_d>0)' ;

zlbshks_f = [n_T(:,zlb_t_f>0)',fg_shock_f] ;
zlbshks_d = [n_T(:,zlb_t_d>0)',fg_shock_d] ;

% Correlations
cc_tpfg_f(:,j_) = vec(corrcoef(fg_shock_save_f(:,j_), tpshock_save_f(:,j_))) ;
cc_tpfg_d(:,j_) = vec(corrcoef(fg_shock_save_d(:,j_), tpshock_save_d(:,j_))) ;

cc_rpfg_f(:,j_) = vec(corrcoef(fg_shock_save_f(:,j_), rpshock_save_f(:,j_))) ;
cc_rpfg_d(:,j_) = vec(corrcoef(fg_shock_save_d(:,j_), rpshock_save_d(:,j_))) ;

cc_zlbshks_f(:,j_) = vec(corrcoef(zlbshks_f)) ;
cc_zlbshks_d(:,j_) = vec(corrcoef(zlbshks_d)) ;

% Calculate term premium shocks
% irf_tp_f      = zeros(n_-1,40) ; 
% irf_tp_f(:,1) = (eye(n_-1)-Pt(:,:,1))\Ct(:,:,1) ; 
% e             = zeros(M_.exo_nbr,40) ; 
% e(1,2)     = -1 ; 
% for t_ = 2 : 40
%    irf_tp_f(:,t_) = Ct(:,:,1) + Pt(:,:,1)*irf_tp_f(:,t_-1) + Dt(:,:,1) * e(:,t_) ; 
% end

end

% === Full set of shock correlations ===
ccloop_save = ccloop ;

ccloop_mean = mean(ccloop,2) ;
ccloop_mean = reshape(ccloop_mean,M_.exo_nbr,M_.exo_nbr) ;

% Remove NANs
idx_nnan    = find(isnan(ccloop_mean(:,1))==0) ; nexo = size(idx_nnan,1) ;
exo_names   = M_.exo_names(idx_nnan,:) ;
ccloop_mean = ccloop_mean(idx_nnan, idx_nnan) ;

idx_nnan_full = find(isnan(ccloop(:,1))==0) ;
ccloop      = ccloop(idx_nnan_full,:) ;

% === Correlations with foreign FG shocks ===
cc_zlbshks_f_mean = mean(cc_zlbshks_f,2) ;
cc_zlbshks_f_mean = reshape(cc_zlbshks_f_mean,M_.exo_nbr + 1,M_.exo_nbr + 1) ; 
idx_nnan    = find(isnan(cc_zlbshks_f_mean(:,1))==0) ; nexo = size(idx_nnan,1) ;
cc_zlbshks_f_mean = cc_zlbshks_f_mean([idx_nnan],[idx_nnan]) ; 

cc_zlbshks_d_mean = mean(cc_zlbshks_d,2) ;
cc_zlbshks_d_mean = reshape(cc_zlbshks_d_mean,M_.exo_nbr + 1,M_.exo_nbr + 1) ; 
idx_nnan    = find(isnan(cc_zlbshks_d_mean(:,1))==0) ; nexo = size(idx_nnan,1) ;
cc_zlbshks_d_mean = cc_zlbshks_d_mean([idx_nnan],[idx_nnan]) ; 

idx_nnan_full = find(isnan(cc_zlbshks_d(:,1))==0) ;
exo_names_zlb = M_.exo_names(idx_nnan(1:end-1),:) ;
cc_zlbshks_f  = cc_zlbshks_f(idx_nnan_full,:) ;
nf            = sqrt(size(cc_zlbshks_f,1)) ; 
cc_zlbshks_f  = cc_zlbshks_f((nf-1)*nf+1 : nf*nf-1,:) ;

cc_zlbshks_d  = cc_zlbshks_d(idx_nnan_full,:) ;
nd            = sqrt(size(cc_zlbshks_d,1)) ; 
cc_zlbshks_d  = cc_zlbshks_d((nd-1)*nd+1 : nd*nd-1,:) ;

%% Plot distributions of correlations

count = 0 ;
for j_ = 1 : nexo  - 1 
    for k_ = j_ + 1 : nexo 

    count = count + 1 ;        
    
    if (count-1)./25 == round((count-1)./25) 
    figure ;
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [20,17*4/5], 'PaperPosition', [0 0 20 17*4/5]) ; 
    set(gcf,'DefaultLineLineWidth', 2) ;
 
    m_ = 1 ;
    end
    
    subplot(5,5,m_) ; 
    [F,Xi] = ksdensity(ccloop((j_-1)*nexo + k_,:)) ;
    plot(Xi,F./sum(F),'b','LineWidth',2) ; hold on ;
    my = ceil(100*max(F./sum(F)))./100 ;
    ylim([0 my]) ;
    yticks([0 : 0.01 : my]) ;
    line([0 0], [0 my],'Color','k','LineWidth',1) ;
    default_graph_options ;    
    
    elmt = strcat([exo_names(j_,:), ', ', exo_names(k_,:)]) ; 
    title(elmt,'Interpreter','none') ;
    
    m_ = m_ + 1 ;
    
    if m_ == 26
       figttle = strcat(['..\figures\ccfig_', num2str(count./25),'.pdf']) ; 
       eval(['print -dpdf ',figttle]) ;
       close
    end
                          
    end
end

print -dpdf '..\figures\ccfig_5' ;

close all ;

%% Plot selected domestic / foreign correlations

foreign_list    = {'eps_xi_f'; 'eps_r_f'; 'eps_g_f'; 'eps_xi_p_f'; 'eps_xi_w_f'; 'eps_r_f_8'} ; 
domestic_list   = {'eps_xi'; 'eps_r'; 'eps_g'; 'eps_xi_H'; 'eps_xi_w'; 'eps_r_8'} ;

ncomp = size(foreign_list,1) ;

figure ;
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [12,10*4/5], 'PaperPosition', [0 0 12 10*4/5]) ; 
set(gcf,'DefaultLineLineWidth', 2) ;

for j_ = 1 : ncomp
fnum = strmatch(foreign_list(j_,:),exo_names,'exact') ;
dnum = strmatch(domestic_list(j_,:),exo_names,'exact') ;

  subplot(2,3,j_) ; 
  [F,Xi] = ksdensity(ccloop((fnum-1)*nexo + dnum,:)) ;
  plot(Xi,F./sum(F),'b','LineWidth',2) ; hold on ;
  my = ceil(100*max(F./sum(F)))./100 ;
  ylim([0 my]) ;
  yticks([0 : 0.01 : my]) ;
  line([0 0], [0 my],'Color','k','LineWidth',1) ;
  default_graph_options ;    
    
  elmt = strcat([exo_names(fnum,:), ', ', exo_names(dnum,:)]) ; 
  title(elmt,'Interpreter','none') ;

end

print -dpdf '..\figures\ccfig_dom_for' ;

close all ;

%% Plot correlations of structural shocks with US FG shocks

figure ;
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [16,14*4/5], 'PaperPosition', [0 0 16 14*4/5]) ; 
set(gcf,'DefaultLineLineWidth', 2) ;

for j_ = 1 : nd  - 1 
      
    subplot(4,4,j_) ; 
    [F,Xi] = ksdensity(cc_zlbshks_f(j_,:)) ;
    plot(Xi,F./sum(F),'b','LineWidth',2) ; hold on ;
    my = ceil(100*max(F./sum(F)))./100 ;
    ylim([0 my]) ;
    yticks([0 : 0.01 : my]) ;
    line([0 0], [0 my],'Color','k','LineWidth',1) ;
    default_graph_options ;    
    
    elmt = strcat([exo_names_zlb(j_,:)]) ; 
    title(elmt,'Interpreter','none') ;
                                
end

print -dpdf '..\figures\cc_FG_f' ;

close all ;

%% Plot correlations of structural shocks with Canadian FG shocks

figure ;
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [20,17*4/5], 'PaperPosition', [0 0 20 17*4/5]) ; 
set(gcf,'DefaultLineLineWidth', 2) ;

for j_ = 1 : nd  - 1 
      
    subplot(4,4,j_) ; 
    [F,Xi] = ksdensity(cc_zlbshks_d(j_,:)) ;
    plot(Xi,F./sum(F),'b','LineWidth',2) ; hold on ;
    my = ceil(100*max(F./sum(F)))./100 ;
    ylim([0 my]) ;
    yticks([0 : 0.01 : my]) ;
    line([0 0], [0 my],'Color','k','LineWidth',1) ;
    default_graph_options ;    
    
    elmt = strcat([exo_names_zlb(j_,:)]) ; 
    title(elmt,'Interpreter','none') ;
                                
end

print -dpdf '..\figures\cc_FG_d' ;

close all ;

%% Correlations of forward guidance shocks and term premia shocks

figure ;
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [12,7*4/5], 'PaperPosition', [0 0 12 7*4/5]) ; 
set(gcf,'DefaultLineLineWidth', 2) ;

subplot(1,2,1) ; 

  [F,Xi] = ksdensity(cc_tpfg_f(2,:)) ;
  plot(Xi,F./sum(F),'b','LineWidth',2) ; hold on ;
  my = ceil(100*max(F./sum(F)))./100 ;
  ylim([0 my]) ;
  yticks([0 : 0.01 : my]) ;
  line([0 0], [0 my],'Color','k','LineWidth',1) ;
  default_graph_options ;    
  title('US FG and TP') ;
  
subplot(1,2,2) ; 

  [F,Xi] = ksdensity(cc_tpfg_d(2,:)) ;
  plot(Xi,F./sum(F),'b','LineWidth',2) ; hold on ;
  my = ceil(100*max(F./sum(F)))./100 ;
  ylim([0 my]) ;
  yticks([0 : 0.01 : my]) ;
  line([0 0], [0 my],'Color','k','LineWidth',1) ;
  default_graph_options ;    
  title('Canada FG and TP') ;

print -dpdf '..\figures\ccfig_FG_TP' ;

close all ;

%% Correlations of forward guidance shocks and consumption risk premium shocks

figure ;
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [12,7*4/5], 'PaperPosition', [0 0 12 7*4/5]) ; 
set(gcf,'DefaultLineLineWidth', 2) ;

subplot(1,2,1) ; 

  [F,Xi] = ksdensity(cc_rpfg_f(2,:)) ;
  plot(Xi,F./sum(F),'b','LineWidth',2) ; hold on ;
  my = ceil(100*max(F./sum(F)))./100 ;
  ylim([0 my]) ;
  yticks([0 : 0.01 : my]) ;
  line([0 0], [0 my],'Color','k','LineWidth',1) ;
  default_graph_options ;    
  title('US FG and xi') ;
  
subplot(1,2,2) ; 

  [F,Xi] = ksdensity(cc_rpfg_d(2,:)) ;
  plot(Xi,F./sum(F),'b','LineWidth',2) ; hold on ;
  my = ceil(100*max(F./sum(F)))./100 ;
  ylim([0 my]) ;
  yticks([0 : 0.01 : my]) ;
  line([0 0], [0 my],'Color','k','LineWidth',1) ;
  default_graph_options ;    
  title('Canada FG and xi') ;

print -dpdf '..\figures\ccfig_FG_RP' ;

close all ;

save simulation_results.mat cc_rpfg_d cc_rpfg_f cc_tpfg_d cc_tpfg_f ccloop_save ;


