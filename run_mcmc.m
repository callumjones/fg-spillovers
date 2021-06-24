%
% run_mcmc.m
% 
% Runs a Markov Chain 
% 
% Jones, Kulish, Rees
% International Spillovers of Forward Guidance Shocks
% Journal of Applied Econometrics, 2021
%

%% Load data

clear all
clc

%%

path(pathdef) ;

warning('Off','all') ;

addpath('./routines')
addpath('./routines/solmat')
addpath('./output') 
addpath('./input')

setup_model

load_data

%% fmincon to get hessian

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

SET.EST.varobs = { ...
    'r_f_1_obs', 'r_f_8_obs', 'dy_f_obs','dc_f_obs', 'infl_f_obs','winf_f_obs', ...
    'r_1_obs', 'r_8_obs', 'dy_obs','dc_obs','infl_obs', 'dy_F_obs', 'dx_obs', 'de_obs', 'winf_obs'} ;

SET.EST.obs = [] ; 
 
for qq=1:length(SET.EST.varobs)
    tmp = find(strcmp(cellstr(SET.M_.endo_names),SET.EST.varobs{qq})) ;
    SET.EST.obs(qq) = tmp ;
end

SET.EST.hidx        = 1:length(SET.EST.obs) ;
SET.EST.hidx_no_r_f = SET.EST.hidx(2:end) ;
SET.EST.hidx_no_r   = SET.EST.hidx([1:6 8:end]) ;
SET.EST.hidx_no_r_f_no_r = SET.EST.hidx([2:6 8:end]) ;

data_start = find(data_dates==1992) ;
zlb_start  = find(data_dates==2008) ;
data       = Zobs(:,data_start:zlb_start-1) ;

SET.nobs = length(data(:,1)) ; 
SET.ss   = length(data(1,:)) ; 

%% minPosterior

x0  = SET.M_.params(SET.EST.params_to_estimate_idx) ;
T_d = zeros(1,length(data)) ;
T_f = zeros(1,length(data)) ;
SET.EST.zlb_at_1pc = zeros(1,SET.ss) ;

P = SET.mats.Q(1:end-1,1:end-1,1) ;
D = SET.mats.G(1:end-1,:,1) ;    

Omega  = eye(SET.variable.l_) ;
DOmega = D(:,:,1)*Omega*D(:,:,1)' ;
SET.Omega = Omega ;
 
Sigma     = (eye((SET.variable.n_-1)^2)-(kron(P,P))) \ DOmega(:) ;
SET.Sigma = reshape(Sigma,SET.variable.n_-1,SET.variable.n_-1) ;  

post_min = @(x) -loglike(x, T_d, T_f, SET, data) - prior(x);

post_min(x0)

% % Commented out and instead load the hessian from previous post min
% [xx_prior,min_param_prior,~,hess1,~,~,~] = ...
%      csminwel(@(x) -prior(x),x0,.005*eye(length(SET.EST.params_to_estimate)),[],1e-6,500);
% 
% [llval_csmw,min_param_ll_csmw,gh,hess_csmw,itct,fcount,retcodeh] = ...
%      csminwel(post_min,x0,hess1,[],1e-5,100);
% 
% hess = inv(hess_csmw) ;
% 
% x0 = min_param_ll_csmw ;
% x_mode_fmin = x0;
% 
% save ./output/hessmat hess x_mode_fmin x0;

load hessmat

%% MCMC

SET.EST.chains = 2 ;
SET.EST.kappa  = 0.2 ;
SET.EST.N      = 100000;
SET.EST.df     = 12 ;
SET.EST.numbersave = 20 ;
SET.aest       = 0 ;
SET.minpost    = 0 ;

data       = Zobs(:,data_start:end) ;
data_dates = data_dates(data_start:end) ;

SET.nobs = length(data(:,1)) ; 
SET.ss   = length(data(1,:)) ;
SET.EST.upperT = 2 ;
SET.EST.Tbstar = 24 ; 

zlb_t_d = zeros(1,SET.ss) ;
zlb_t_f = zeros(1,SET.ss) ;

% Canadian effective percent lower bound ends on 2015Q1
zlb_t_d(1,find(data_dates==2009.25):find(data_dates==2014.75)) = 1 ;

% US effective lower bound ends on 2015Q4
zlb_t_f(1,find(data_dates==2009):find(data_dates==2015.5))     = 1 ;
zlb_t = [zlb_t_d ; zlb_t_f] ;

SET.EST.zlb_at_1pc = zeros(1,SET.ss) ;
SET.EST.zlb_at_1pc(1,find(data_dates==2010.5):find(data_dates==2014.75))=1 ;

c = SET.EST.kappa * inv(hess) ; % Use fmincon Hessian
c = chol(c)' ;

maxproc = SET.EST.chains ;
%parpool(maxproc)

clear params_x params_T_d params_T_f  acc_rates reject_rate
clc

parfor runs=1:maxproc
    warning('Off','all') ;

    out = mcmc(x0, zlb_t, c, SET, data, runs) ;
    
    params_x(:,:,runs)   = out.xx_out ;
    params_T_d(:,:,runs) = out.Tbs_out_d ;
    params_T_f(:,:,runs) = out.Tbs_out_f ;
    acc_rates(runs,:)    = out.acc_rates ;

end

clear date

eval(['save ./output/mhall_', date, '_estimatechis.mat;']);

%% Stack the (trimmed) chain outputs

params_x_   = [] ;
params_T_d_ = [] ;
params_T_f_ = [] ;

burn = 50000 ;

for runs=1:maxproc
    params_x_ = [params_x_ ; params_x(burn:end,:,runs)] ; 
    params_T_d_ = [params_T_d_ ; params_T_d(burn:end,:,runs)] ;
    params_T_f_ = [params_T_f_ ; params_T_f(burn:end,:,runs)] ; 
end

%% Chain diagnostics

sq_R_x = [];
sq_R_T_d = [];
sq_R_T_f = [];

pvec=1:52;

count = 1 ;

yvec = [200 300 400 500 600 700 800 900 1000 1250 1500 2000 3000 4000] ;
yvec = [yvec 5000:10000:length(params_x(:,1,1))] ;

for j=yvec
    j
    sq_R_x(:,count) = ...
        mh_diag(params_x(1:j,pvec,:)) ;
    if ~isempty(params_T_d)
    sq_R_T_d(:,count) = ...
        mh_diag(params_T_d(1:j,:,:)) ;
    end
    if ~isempty(params_T_f)
    sq_R_T_f(:,count) = ...
        mh_diag(params_T_f(1:j,:,:)) ;
    end
    count = count+1 ;
end

figure ;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 8]) ;
subplot(3,1,1) ; hold on; 
    plot(log(yvec),sq_R_x','k','LineWidth',1) ;
    line([log(yvec(1)) log(yvec(end))],[1.1 1.1],'LineWidth',1,'color','r','LineStyle','--') ;
    line([log(yvec(1)) log(yvec(end))],[1 1],'LineWidth',1,'color','k') ;
    box('on') ; grid ;
    title('Parameters') ;
    xlim([log(yvec(1)) log(yvec(end))]) ;

subplot(3,1,2) ; hold on; 
    plot(log(yvec),sq_R_T_d','k','LineWidth',1) ;
    line([log(yvec(1)) log(yvec(end))],[1.1 1.1],'LineWidth',1,'color','r','LineStyle','--') ;
    line([log(yvec(1)) log(yvec(end))],[1 1],'LineWidth',1,'color','k') ;
    box('on') ; grid ;
    ylabel('R^2 diagnostic') ;
    title('Canadian Durations') ;
    xlim([log(yvec(1)) log(yvec(end))]) ;

subplot(3,1,3) ; hold on; 
    plot(log(yvec),sq_R_T_f','k','LineWidth',1) ;
    line([log(yvec(1)) log(yvec(end))],[1.1 1.1],'LineWidth',1,'color','r','LineStyle','--') ;
    line([log(yvec(1)) log(yvec(end))],[1 1],'LineWidth',1,'color','k') ;
    box('on') ; grid ;
    xlabel('Length of chain, log scale') ;
    title('US Durations') ;
    xlim([log(yvec(1)) log(yvec(end))]) ;
set(findall(gcf,'-property','FontSize'),'FontSize',9) ;
print -depsc ./output/gelman_diagnostics.eps
close

%% Plot Data

if 0
    
figure ;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 10]) ;
set(gcf,'DefaultLineLineWidth', 2) ;
subplot(4,4,1) ;
    plot(data_dates,1*(data(3,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
   % ylim([-3 2]) ;
    xlim([1992 2020]) ;
    title('US Output, \% $\Delta$') ;   
subplot(4,4,2) ;
    plot(data_dates,1*(data(4,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
   % ylim([0 8]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    xlim([1992 2020]) ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
	title('US Consumption, \% $\Delta$') ;
subplot(4,4,3) ;
    plot(data_dates,1*(data(1,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    %ylim([-2 15]) ;
    ax = gca ;
    ax.YGrid='on';  
    xlim([1992 2020]) ;
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
    title('Fed Funds Rate, \%') ;    
subplot(4,4,4) ;
    plot(data_dates,1*(data(5,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    %ylim([-2 15]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    xlim([1992 2020]) ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
    title('US Inflation, \%') ;        
subplot(4,4,5) ;
    plot(data_dates,1*(data(6,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    %ylim([0 8]) ;
    ax = gca ;
    xlim([1992 2020]) ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
	title('US Nominal Wage \% $\Delta$') ;    
subplot(4,4,6) ;
    plot(data_dates,1*(data(9,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    xlim([1992 2020]) ;
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
   % ylim([-3 2]) ;
    title('Canada Output, \% $\Delta$') ;   
subplot(4,4,7) ;
    plot(data_dates,1*(data(10,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    xlim([1992 2020]) ;
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
   % ylim([0 8]) ;
	title('Canada Consumption, \% $\Delta$') ; 
subplot(4,4,8) ;
    plot(data_dates,1*(data(7,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    %ylim([-2 15]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    xlim([1992 2020]) ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
    title('Canada Bank Rate, \%') ;    
subplot(4,4,9) ;
    plot(data_dates,1*(data(11,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    %ylim([-2 15]) ;
    ax = gca ;
    ax.YGrid='on';
    xlim([1992 2020]) ;
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
    title('Canada Inflation, \%') ;      
subplot(4,4,10) ;
    plot(data_dates,1*data(15,:),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    xlim([1992 2020]) ;
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
	title('Canada Nominal Wage, \% $\Delta$') ;     
subplot(4,4,11) ;
    plot(data_dates,1*data(12,:),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    xlim([1992 2020]) ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
	title('Imports, \% $\Delta$') ; 
subplot(4,4,12) ;
    plot(data_dates,1*data(13,:),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    xlim([1992 2020]) ;    
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
	title('Exports, \% $\Delta$') ;     
subplot(4,4,13) ;
    plot(data_dates,1*data(14,:),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    xlim([1992 2020]) ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
	title('USD/CAD Exchange Rate, \% $\Delta$') ;         
subplot(4,4,14) ;
    plot(data_dates,1*(data(2,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    %ylim([-2 15]) ;
    ax = gca ;
    ax.YGrid='on';
    xlim([1992 2020]) ;
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
	title('US 2Y Interest Rate, \%') ; 
subplot(4,4,15) ;
    plot(data_dates,1*(data(8,:)),'k','LineWidth',2)
    line([1980 2020],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    xlim([1992 2020]) ;
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
    %ylim([-2 15]) ;
	title('Canada 2Y Interest Rate, \%') ;   
set(findall(gcf,'-property','FontSize'),'FontSize',11) ;
print -depsc ./output/data.eps
close ;   

end