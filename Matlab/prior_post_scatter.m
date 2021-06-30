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

setup_model

%%

load mhall_05-Mar-2020_estimatechis.mat

num_params = length(SET.EST.params_to_estimate) ;

%priorx = draw_from_prior(SET, params_x(end,:,1)) ;
load prior_sav 

params_x_   = [] ;
params_T_d_ = [] ;
params_T_f_ = [] ;

burn = 50000 ;

for runs=1:maxproc
    params_x_ = [params_x_ ; params_x(burn:end,:,runs)] ; 
    params_T_d_ = [params_T_d_ ; params_T_d(burn:end,:,runs)] ;
    params_T_f_ = [params_T_f_ ; params_T_f(burn:end,:,runs)] ; 
end

%%

elas_vec = linspace(eps,5,20) ;

T = 30 ;

et = 0*randn(SET.variable.l_,T) ;
et(SET.variable.shock.eps_r_f,2) = -1 ;

NN=1500 ;

rand_prior = randsample(1:length(priorx),NN) ;
rand_post  = randsample(1:length(params_x_),NN) ;

%% Prior

%rand_post = randsample(1:length(priorx),NN) ;

counter = 1 ;

for jjj=1:NN
    
    jjj

    dyn_in.M_.params(SET.EST.params_to_estimate_idx) = priorx(rand_prior(jjj),:) ;
    dursest_f = zeros(1,T) ;
    dursest_d = zeros(1,T) ;

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
    
    if isempty(SET.mats.Q) ; continue ; end

    SET.y_ss = (eye(n_-1) - SET.mats.Q(1:n_-1,1:n_-1)) \ SET.mats.Q(1:n_-1,end)  ;
    SET.y_ss = [SET.y_ss ; 1] ; % Constant
    
    y_zlb = SET.y_ss ;
    z_f_T = zeros(1,T) ; 
    z_d_T = zeros(1,T) ; 

    for t = 2:T
        [y_zlb(:,t), ~, ~, z_f_T(t), z_d_T(t)] = ...
            zlb(SET, y_zlb(:,t-1), et(:,t)) ;        
    end

    irfs_prior(:,:,counter) = y_zlb ;

    elast_vec_pr(counter) = dyn_in.M_.params(SET.param_names.tau) ;
    
    counter = counter+1 ;

end

%% Posterior

counter = 1 ;

for jjj=1:NN
    
    jjj

    dyn_in.M_.params(SET.EST.params_to_estimate_idx) = params_x_(rand_post(jjj),:) ;
    dursest_f = zeros(1,T) ;
    dursest_d = zeros(1,T) ;

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
    
    if isempty(SET.mats.Q) ; continue ; end

    SET.y_ss = (eye(n_-1) - SET.mats.Q(1:n_-1,1:n_-1)) \ SET.mats.Q(1:n_-1,end)  ;
    SET.y_ss = [SET.y_ss ; 1] ; % Constant
    
    y_zlb = SET.y_ss ;
    z_f_T = zeros(1,T) ; 
    z_d_T = zeros(1,T) ; 

    for t = 2:T
        [y_zlb(:,t), ~, ~, z_f_T(t), z_d_T(t)] = ...
            zlb(SET, y_zlb(:,t-1), et(:,t)) ;        
    end

    irfs_post(:,:,counter) = y_zlb ;

    elast_vec_post(counter) = dyn_in.M_.params(SET.param_names.tau) ;
    
    counter = counter+1 ;

end

return

%%

figure ;
set(gcf,'DefaultLineLineWidth', 3) ;
hold on ;
    plot(elast_vec_pr(1),100*irfs_prior(SET.variable.y,2,1),'-o','Color',[0 127/256 190/256]);
    plot(elast_vec_post(1),100*irfs_post(SET.variable.y,2,1),'-ok');
    
    scatter(elast_vec_pr(:),squeeze(100*irfs_prior(SET.variable.y,2,:)),150,'.','MarkerFaceColor',[0 127/256 190/256],'MarkerEdgeColor',[0 127/256 190/256]);
    scatter(elast_vec_post(:),squeeze(100*irfs_post(SET.variable.y,2,:)),150,'.k');
    %title('Domestic Output Response','Interpreter','latex')
    line([0 4],[0 0],'LineWidth',1,'Color','k') ;
    ylim([-0.15 0.0500001]) ;
    xlim([0 4]) ;
    %ylim([-0.3 0.2]) ;
    %box('on') ; 
    %grid('on') ;    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',20,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')     
        
    xlabel('$\tau$','Interpreter','latex') ;
    h=legend('Draws From Prior','Draws From Posterior','location','northeast') ;
    set(h,'Interpreter','latex','FontSize',20);  
%set(findall(gcf,'-property','FontSize'),'FontSize',11) ;
print -depsc ./output/fig_prior_post_over_tau.eps
close ; 

%%

for ii=1:length(params_x_(1,:))
    %[N,X] = hist(params_x_(:,ii),5000) ;
    [a,b]=ksdensity(params_x_(:,ii),'NumPoints',200);
    %tmp = X(N==max(N)) ;
    x_mode_fmin(ii) = b(a==max(a)) ;
end

dyn_in.M_.params(SET.EST.params_to_estimate_idx) = x_mode_fmin ;

dyn_in.M_.params(SET.param_names.tau) = 0.5 ;

dursest_f = zeros(1,T) ;
dursest_d = zeros(1,T) ;

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

SET.y_ss = (eye(n_-1) - SET.mats.Q(1:n_-1,1:n_-1)) \ SET.mats.Q(1:n_-1,end)  ;
SET.y_ss = [SET.y_ss ; 1] ; % Constant

y_zlb = SET.y_ss ;
z_f_T = zeros(1,T) ; 
z_d_T = zeros(1,T) ; 

for t = 2:T
    [y_zlb(:,t), ~, ~, z_f_T(t), z_d_T(t)] = ...
        zlb(SET, y_zlb(:,t-1), et(:,t)) ;        
end
    
irfs_lowtau = y_zlb ;

dyn_in.M_.params(SET.param_names.tau) = 5 ;

dursest_f = zeros(1,T) ;
dursest_d = zeros(1,T) ;

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

SET.y_ss = (eye(n_-1) - SET.mats.Q(1:n_-1,1:n_-1)) \ SET.mats.Q(1:n_-1,end)  ;
SET.y_ss = [SET.y_ss ; 1] ; % Constant

y_zlb = SET.y_ss ;
z_f_T = zeros(1,T) ; 
z_d_T = zeros(1,T) ; 

for t = 2:T
    [y_zlb(:,t), ~, ~, z_f_T(t), z_d_T(t)] = ...
        zlb(SET, y_zlb(:,t-1), et(:,t)) ;        
end

irfs_hightau = y_zlb ;

%%


figure; 
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 5]) ;
set(gcf,'DefaultLineLineWidth', 2) ;
subplot(2,2,1) ; hold on;
    plot(100*(irfs_lowtau(SET.variable.y_f,:)),'k','LineWidth',2);  
    plot(100*(irfs_hightau(SET.variable.y_f,:)),'-.','LineWidth',2); 
    line([0 30],[0 0], 'Color', 'k', 'LineWidth',1) ;
    %plot(tmp,100*(exp(x_t(SET.variable.iobs,:)).^4-1)) ;
    %xlim([1990 2020]) ;

    %ylim([0 10]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 

    h = legend('Low $\tau$', 'High $\tau$', 'Location','NorthEast');
    set(h,'Interpreter','latex','FontSize',10);      
    
    title('A. US Output, \%') ;
subplot(2,2,2) ; hold on;
    plot(1*(irfs_lowtau(SET.variable.r_f_1_obs,:)),'k','LineWidth',2);  
    plot((irfs_hightau(SET.variable.r_f_1_obs,:)),'-.','LineWidth',2);  
    %plot(tmp,100*(exp(x_t(SET.variable.iobs,:)).^4-1)) ;
    %xlim([1990 2020]) ;

    %ylim([0 10]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 

    %h = legend('Low $\tau$', 'High $\tau$', 'Location','NorthEast');
    %set(h,'Interpreter','latex','FontSize',10);      
    
    title('B. Fed Funds Rate, \%') ;  
subplot(2,2,3) ; hold on;
line([0 30],[0 0], 'Color', 'k', 'LineWidth',1) ;
    plot(100*(irfs_lowtau(SET.variable.y,:)),'k','LineWidth',2);  
    plot(100*(irfs_hightau(SET.variable.y,:)),'-.','LineWidth',2);  
    %plot(tmp,100*(exp(x_t(SET.variable.iobs,:)).^4-1)) ;
    %xlim([1990 2020]) ;

    %ylim([0 10]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 

    %h = legend('Low $\tau$', 'High $\tau$', 'Location','NorthEast');
    %set(h,'Interpreter','latex','FontSize',10);      
    
    title('C. Canada Output, \%') ;
subplot(2,2,4) ;hold on;
    
    plot((irfs_lowtau(SET.variable.r_1_obs,:)),'k','LineWidth',2);  
    plot((irfs_hightau(SET.variable.r_1_obs,:)),'-.','LineWidth',2);  
    %plot(tmp,100*(exp(x_t(SET.variable.iobs,:)).^4-1)) ;
    %xlim([1990 2020]) ;

    %ylim([0 10]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',11,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
    
    title('D. Canada Bank Rate, \%') ;  
set(findall(gcf, 'Type', 'Text'),'FontWeight', 'Normal')
%set(h,'FontSize',10) ;
print -depsc ./output/irf_lowhigh_tau.eps
close  


