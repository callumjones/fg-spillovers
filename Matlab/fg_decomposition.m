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

load_data

%%

load mhall_05-Mar-2020_estimatechis.mat

SET.maxchk   = 50 ; % Check ZLB binding out 50 periods

% Select number of posterior draws to do counterfactuals for
% To generate the figures in the paper, use max_id = 10000
max_id = 100 ; 

params_x_   = [] ;
params_T_d_ = [] ;
params_T_f_ = [] ;

burn = 50000 ;

for runs=1:maxproc
    params_x_ = [params_x_ ; params_x(burn:end,:,runs)] ; 
    params_T_d_ = [params_T_d_ ; params_T_d(burn:end,:,runs)] ;
    params_T_f_ = [params_T_f_ ; params_T_f(burn:end,:,runs)] ; 
end

ind = randsample(1:length(params_x_),max_id) ;

eps_xi_f_vec=[];
lb_T_vec=[];
lb_T_d_vec=[];
z_f_T_fg_zlb_vec=[];
z_d_T_fg_zlb_vec=[];

%SET.EST.zlb_at_1pc = 0*SET.EST.zlb_at_1pc ;

for jj=1:max_id
    
ii=ind(jj) ;
    
dyn_in.M_.params(SET.EST.params_to_estimate_idx) = params_x_(ii,:) ;
dursest_f = params_T_f_(ii,:) ;
dursest_d = params_T_d_(ii,:) ;

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

SET.nobs = length(data(:,1)) ; 
SET.ss   = length(data(1,:)) ; 

ss      = SET.ss ;
n       = SET.variable.n_ ; 
nparam  = SET.nparam ;
nobs    = SET.nobs ;

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

SET.y_ss_init = (eye(n-1) - Qt(1:n-1,1:n-1,1)) \ Qt(1:n-1,end,1)  ;
SET.y_ss_init = [SET.y_ss_init ; 1] ;

%%% Matrices for US FG Only
SC.t_f = zeros(1,SET.ss) ;
SC.t_d = zeros(1,SET.ss) ;  
SC.t_f(zlb_t_f>0) = dursest_f ; tmp_FG_US = SC.t_f ;
mats_fg_US = tv_mats(SET,SC) ;

Qt_fg_US = mats_fg_US.Qt ;
Gt_fg_US = mats_fg_US.Gt ;
%%%

%%% Matrices for Can FG Only
SC.t_f = zeros(1,SET.ss) ;
SC.t_d = zeros(1,SET.ss) ;  
SC.t_d(zlb_t_d>0) = dursest_d ; tmp_FG_Can = SC.t_d ;
mats_fg_Can = tv_mats(SET,SC) ;

Qt_fg_Can = mats_fg_Can.Qt ;
Gt_fg_Can = mats_fg_Can.Gt ;
%%%

y_no_sto  = SET.y_ss_init ;
y_no_zlb  = x_T(:,1) ;
y_zlb     = x_T(:,1) ;
y_fg_US   = x_T(:,1) ;
y_fg_Can  = x_T(:,1) ;
z_f_T     = zeros(1,SET.ss) ;
z_d_T     = zeros(1,SET.ss) ;
z_f_T_fg_US = zeros(1,SET.ss) ;
z_d_T_fg_US = zeros(1,SET.ss) ;
z_f_T_fg_Can = zeros(1,SET.ss) ;
z_d_T_fg_Can = zeros(1,SET.ss) ;

Qt_nozlb = repmat(Qt(:,:,1),[1 1 SET.ss+SET.maxchk]) ;
Gt_nozlb = repmat(Gt(:,:,1),[1 1 SET.ss+SET.maxchk]) ;

Qt_fg_US(:,:,(SET.ss+1):(SET.ss+SET.maxchk)) = repmat(Qt(:,:,1),[1 1 SET.maxchk]) ;
Gt_fg_US(:,:,(SET.ss+1):(SET.ss+SET.maxchk)) = repmat(Gt(:,:,1),[1 1 SET.maxchk]) ;

Qt_fg_Can(:,:,(SET.ss+1):(SET.ss+SET.maxchk)) = repmat(Qt(:,:,1),[1 1 SET.maxchk]) ;
Gt_fg_Can(:,:,(SET.ss+1):(SET.ss+SET.maxchk)) = repmat(Gt(:,:,1),[1 1 SET.maxchk]) ;

ychk = x_T(:,1) ;

for t=2:SET.ss
    y_no_sto(:,t) = SET.mats.Q*y_no_sto(:,t-1) ;
    y_no_zlb(:,t) = SET.mats.Q*y_no_zlb(:,t-1) + SET.mats.G*n_T(:,t) ;
    
    ychk(:,t) = Qt(:,:,t)*ychk(:,t-1) + Gt(:,:,t)*n_T(:,t) ;
    
    [y_zlb(:,t), ~, ~, z_f_T(t), z_d_T(t)] = ...
        zlb(SET, y_zlb(:,t-1), n_T(:,t), Qt_nozlb, Gt_nozlb, t) ;
   
    FG.z_f_i = (tmp_FG_US(t)>0) + 100*(tmp_FG_US(t)==0) ; 
    FG.z_f_f = tmp_FG_US(t) ; 
    FG.z_d_i = 100 ; 
    FG.z_d_f = 0 ; 

    [y_fg_US(:,t), ~, ~, z_f_T_fg_US(t), z_d_T_fg_US(t)] = ...
        zlb_usfg(SET, y_fg_US(:,t-1), n_T(:,t), t, FG) ;
    
    [y_fg_Can(:,t), ~, ~, z_f_T_fg_Can(t), z_d_T_fg_Can(t)] = ...
        zlb(SET, y_fg_Can(:,t-1), n_T(:,t), Qt_fg_Can, Gt_fg_Can, t) ;
end

eps_xi_f_vec(jj,:) = n_T(SET.variable.shock.eps_xi_f,:) ;
eps_r_f_vec(jj,:)  = n_T(SET.variable.shock.eps_r_f,:) ;
eps_r_vec(jj,:)    = n_T(SET.variable.shock.eps_r,:) ;
lb_T_vec(jj,:)     = z_f_T(zlb_t_f>0) ; % ZLB duration
lb_T_d_vec(jj,:)   = z_d_T(zlb_t_d>0) ; % ZLB duration
y_zlb_mat(:,:,jj)  = y_zlb([SET.variable.r_f_1_obs SET.variable.r_1_obs SET.variable.infl_f_obs SET.variable.infl_obs],:) ;
y_fg_US_mat(:,:,jj)  = y_fg_US([SET.variable.r_f_1_obs SET.variable.r_1_obs SET.variable.infl_f_obs SET.variable.infl_obs],:) ;
y_fg_Can_mat(:,:,jj) = y_fg_Can([SET.variable.r_f_1_obs SET.variable.r_1_obs SET.variable.infl_f_obs SET.variable.infl_obs],:) ;
y_no_zlb_mat(:,:,jj) = y_no_zlb([SET.variable.r_f_1_obs SET.variable.r_1_obs SET.variable.infl_f_obs SET.variable.infl_obs],:) ;

output_f_no_fg(:,jj) = gen_idx(y_zlb(SET.variable.dy_f_obs,:)/100) ;
output_no_fg(:,jj)   = gen_idx(y_zlb(SET.variable.dy_obs,:)/100) ;

output_f_fg_US(:,jj) = gen_idx(y_fg_US(SET.variable.dy_f_obs,:)/100) ;
output_fg_US(:,jj)   = gen_idx(y_fg_US(SET.variable.dy_obs,:)/100) ;

output_f_fg_Can(:,jj) = gen_idx(y_fg_Can(SET.variable.dy_f_obs,:)/100) ;
output_fg_Can(:,jj)   = gen_idx(y_fg_Can(SET.variable.dy_obs,:)/100) ;

output_f_nozlb(:,jj) = gen_idx(y_no_zlb(SET.variable.dy_f_obs,:)/100) ;
output_nozlb(:,jj)   = gen_idx(y_no_zlb(SET.variable.dy_obs,:)/100) ;

output_f       = gen_idx(x_T(SET.variable.dy_f_obs,:)/100) ;
output         = gen_idx(x_T(SET.variable.dy_obs,:)/100) ;

output_f_chk(:,jj) = gen_idx(ychk(SET.variable.dy_f_obs,:)/100) ;
output_chk(:,jj) = gen_idx(ychk(SET.variable.dy_obs,:)/100) ;

fprintf('\n')
fprintf('\n')
fprintf('Draw %i \n', jj)
fprintf('\n')
fprintf('  Cumulative output loss if no FG:  US     = %9.2f \n', sum(100*((output_f_no_fg(2:end,jj)./output_f(2:end)')-1)));
fprintf('  Cumulative output loss if no FG:  Canada = %9.2f \n', sum(100*((output_no_fg(2:end-1,jj)./output(2:end-1)')-1)));

cumulative_y_loss_us(jj) = sum(100*((output_f_no_fg(2:end,jj)./output_f(2:end)')-1)) ;
cumulative_y_loss_can(jj) = sum(100*((output_no_fg(2:end-1,jj)./output(2:end-1)')-1)) ;

end

%save ./output/decomp_sav

%% Draw Counterfactual Figures

% Randomly pick max_id_plot plots to draw
max_id_plot = max_id ;
ind_plot = randsample(1:max_id,max_id_plot) ;

SwatheOpt = PlotSwatheOption;
SwatheOpt.marker = '*';
SwatheOpt.trans = 1;
SwatheOpt.xaxis =data_dates';
%SwatheOpt.linecol = [0.8500, 0.3250, 0.0980];

%% No FG in both US and Canada

figure ; 

set(gcf, 'DefaultAxesLineWidth', 1);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
set(gcf, 'DefaultLegendInterpreter','latex');

subplot(2,2,1) ; hold on ;
    plot(data_dates,4*x_T(SET.variable.infl_f_obs,:),'LineWidth',2,'Color','k') ;
    PlotSwathe(4*squeeze(mean(y_zlb_mat(3,:,ind_plot),3)),...
        4*squeeze(prctile(y_zlb_mat(3,:,ind_plot),[5 95],3)),SwatheOpt); hold on;
	plot(data_dates,4*x_T(SET.variable.infl_f_obs,:),'LineWidth',2,'Color','k') ;

    xlim([2006 2019]) ;
    title('A. US, \% Annual','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
  
subplot(2,2,2) ; hold on ;
    plot(data_dates,4*x_T(SET.variable.infl_obs,:),'LineWidth',2,'Color','k') ;
    PlotSwathe(4*squeeze(mean(y_zlb_mat(4,:,ind_plot),3)),...
        4*squeeze(prctile(y_zlb_mat(4,:,ind_plot),[5 95],3)),SwatheOpt); hold on;
    plot(data_dates,4*x_T(SET.variable.infl_obs,:),'LineWidth',2,'Color','k') ;

    line([2006 2019],[0 0],'LineWidth',1,'Color','k') ;
    
    xlim([2006 2019]) ;
    title('B. Canada, \% Annual','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;

print -depsc ./output/cf_nofg_infl_chipchiw_v4.eps
close  


figure ; 

set(gcf, 'DefaultAxesLineWidth', 1);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
set(gcf, 'DefaultLegendInterpreter','latex');

subplot(2,2,1) ; hold on ;
    plot(data_dates,mean(output_f_chk'),'LineWidth',2,'Color','k') ;
    PlotSwathe(squeeze(mean(output_f_no_fg(:,ind_plot),2)),squeeze(prctile(output_f_no_fg(:,ind_plot),[5 95],2)),SwatheOpt); hold on;

    plot(data_dates,mean(output_f_chk'),'LineWidth',2,'Color','k') ;
    xlim([2006 2019]) ;
    title('A. US, Index (1992=1)','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    h=legend('Data','Counterfactuals','location','northwest');
	set(h,'Interpreter','latex','FontSize',8);    
    
subplot(2,2,2) ; hold on ;
    plot(data_dates,mean(output_chk'),'LineWidth',2,'Color','k') ;
    PlotSwathe(squeeze(mean(output_no_fg(:,ind_plot),2)),squeeze(prctile(output_no_fg(:,ind_plot),[5 95],2)),SwatheOpt); hold on;
    plot(data_dates,mean(output_chk'),'LineWidth',2,'Color','k') ;

    xlim([2006 2019]) ;
    title('B. Canada, Index (1992=1)','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;

print -depsc ./output/cf_nofg_y_chipchiw_v4.eps
close  

tmp=100*(log(mean(output_f_chk'))-log(squeeze(mean(output_f_no_fg(:,ind_plot),2)')));
mean(tmp(find(data_dates==2009):find(data_dates==2015.5)))
tmp=100*(log(mean(output_chk'))-log(squeeze(mean(output_no_fg(:,ind_plot),2)')));
mean(tmp(find(data_dates==2009):find(data_dates==2015.5)))

%% FG in Canada Only

figure ; 

set(gcf, 'DefaultAxesLineWidth', 1);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
set(gcf, 'DefaultLegendInterpreter','latex');

subplot(2,2,1) ; hold on ;
    plot(data_dates,4*x_T(SET.variable.infl_f_obs,:),'LineWidth',2,'Color','k') ;
    PlotSwathe(4*squeeze(mean(y_fg_Can_mat(3,:,ind_plot),3)),...
        4*squeeze(prctile(y_fg_Can_mat(3,:,ind_plot),[5 95],3)),SwatheOpt); hold on;

    plot(data_dates,4*x_T(SET.variable.infl_f_obs,:),'LineWidth',2,'Color','k') ;

    xlim([2006 2019]) ;
    title('A. US, \% Annual','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    
subplot(2,2,2) ; hold on ;
    plot(data_dates,4*x_T(SET.variable.infl_obs,:),'LineWidth',2,'Color','k') ;

    PlotSwathe(4*squeeze(mean(y_fg_Can_mat(4,:,ind_plot),3)),...
        4*squeeze(prctile(y_fg_Can_mat(4,:,ind_plot),[5 95],3)),SwatheOpt); hold on;

    plot(data_dates,4*x_T(SET.variable.infl_obs,:),'LineWidth',2,'Color','k') ;

    line([2006 2019],[0 0],'LineWidth',1,'Color','k') ;
    
    xlim([2006 2019]) ;
    ylim([-1 4]) ;
    title('B. Canada, \% Annual','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;

print -depsc ./output/cf_fg_Can_infl_chipchiw_v4.eps
close  


figure ; 

set(gcf, 'DefaultAxesLineWidth', 1);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
set(gcf, 'DefaultLegendInterpreter','latex');

subplot(2,2,1) ; hold on ;
    plot(data_dates,mean(output_f_chk'),'LineWidth',2,'Color','k') ;
    PlotSwathe(squeeze(mean(output_f_fg_Can(:,ind_plot),2)),squeeze(prctile(output_f_fg_Can(:,ind_plot),[5 95],2)),SwatheOpt); hold on;

    plot(data_dates,mean(output_f_chk'),'LineWidth',2,'Color','k') ;
    xlim([2006 2019]) ;
    title('A. US, Index (1992=1)','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    h=legend('Data','Counterfactuals','location','northwest');
	set(h,'Interpreter','latex','FontSize',8);    
    
subplot(2,2,2) ; hold on ;
    plot(data_dates,mean(output_chk'),'LineWidth',2,'Color','k') ;
    PlotSwathe(squeeze(mean(output_fg_Can(:,ind_plot),2)),squeeze(prctile(output_fg_Can(:,ind_plot),[5 95],2)),SwatheOpt); hold on;
    plot(data_dates,mean(output_chk'),'LineWidth',2,'Color','k') ;

    xlim([2006 2019]) ;
    title('B. Canada, Index (1992=1)','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;

print -depsc ./output/cf_fg_Can_y_chipchiw_v4.eps
close  

tmp=100*(log(mean(output_chk'))-log(squeeze(mean(output_fg_Can(:,ind_plot),2)')));
mean(tmp(find(data_dates==2009):find(data_dates==2015.5)))

tmp= 4*x_T(SET.variable.infl_obs,:)-4*squeeze(mean(y_fg_Can_mat(4,:,ind_plot),3));
mean(tmp(find(data_dates==2009):find(data_dates==2015.5)))

%% FG in US Only

figure ; 

set(gcf, 'DefaultAxesLineWidth', 1);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
set(gcf, 'DefaultLegendInterpreter','latex');

subplot(2,2,1) ; hold on ;
    plot(data_dates,4*x_T(SET.variable.infl_f_obs,:),'LineWidth',2,'Color','k') ;
    PlotSwathe(4*squeeze(mean(y_fg_US_mat(3,:,ind_plot),3)),...
        4*squeeze(prctile(y_fg_US_mat(3,:,ind_plot),[5 95],3)),SwatheOpt); hold on;
	plot(data_dates,4*x_T(SET.variable.infl_f_obs,:),'LineWidth',2,'Color','k') ;

    xlim([2006 2019]) ;
    title('A. US, \% Annual','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    
subplot(2,2,2) ; hold on ;
    plot(data_dates,4*x_T(SET.variable.infl_obs,:),'LineWidth',2,'Color','k') ;
    PlotSwathe(4*squeeze(mean(y_fg_US_mat(4,:,ind_plot),3)),...
        4*squeeze(prctile(y_fg_US_mat(4,:,ind_plot),[5 95],3)),SwatheOpt); hold on;

    plot(data_dates,4*x_T(SET.variable.infl_obs,:),'LineWidth',2,'Color','k') ;

    line([2006 2019],[0 0],'LineWidth',1,'Color','k') ;
    
    xlim([2006 2019]) ;
    ylim([-1 4]) ;
    title('B. Canada, \% Annual','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;

print -depsc ./output/cf_fg_US_infl_chipchiw_v4.eps
close  


figure ; 

set(gcf, 'DefaultAxesLineWidth', 1);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
set(gcf, 'DefaultLegendInterpreter','latex');

subplot(2,2,1) ; hold on ;
    plot(data_dates,mean(output_f_chk'),'LineWidth',2,'Color','k') ;
    PlotSwathe(squeeze(mean(output_f_fg_US(:,ind_plot),2)),squeeze(prctile(output_f_fg_US(:,ind_plot),[5 95],2)),SwatheOpt); hold on;

    plot(data_dates,mean(output_f_chk'),'LineWidth',2,'Color','k') ;
    xlim([2006 2019]) ;
    title('A. US, Index (1992=1)','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    h=legend('Data','Counterfactuals','location','northwest');
    %legend('boxoff')
	set(h,'Interpreter','latex','FontSize',8);    
    
subplot(2,2,2) ; hold on ;
    plot(data_dates,mean(output_chk'),'LineWidth',2,'Color','k') ;
    PlotSwathe(squeeze(mean(output_fg_US(:,ind_plot),2)),squeeze(prctile(output_fg_US(:,ind_plot),[5 95],2)),SwatheOpt); hold on;
    plot(data_dates,mean(output_chk'),'LineWidth',2,'Color','k') ;

    xlim([2006 2019]) ;
    title('B. Canada, Index (1992=1)','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;

print -depsc ./output/cf_fg_US_y_chipchiw_v4.eps
close  


figure ; 

set(gcf, 'DefaultAxesLineWidth', 1);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
set(gcf, 'DefaultLegendInterpreter','latex');

subplot(2,2,1) ; hold on ;
    plot(data_dates,mean(output_chk'),'LineWidth',2,'Color','k') ;
    PlotSwathe(squeeze(mean(output_fg_US(:,ind_plot),2)),squeeze(prctile(output_fg_US(:,ind_plot),[5 95],2)),SwatheOpt); hold on;
    plot(data_dates,mean(output_chk'),'LineWidth',2,'Color','k') ;

    xlim([2006 2019]) ;
    title('A. Output, Index (1992=1)','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	h=legend('Data','Counterfactuals','location','northwest');
	set(h,'Interpreter','latex','FontSize',7);    

subplot(2,2,2) ; hold on ;
    plot(data_dates,4*x_T(SET.variable.infl_obs,:),'LineWidth',2,'Color','k') ;
    PlotSwathe(4*squeeze(mean(y_fg_US_mat(4,:,ind_plot),3)),...
        4*squeeze(prctile(y_fg_US_mat(4,:,ind_plot),[5 95],3)),SwatheOpt); hold on;

    plot(data_dates,4*x_T(SET.variable.infl_obs,:),'LineWidth',2,'Color','k') ;

    line([2006 2019],[0 0],'LineWidth',1,'Color','k') ;
    
    xlim([2006 2019]) ;
    ylim([-1 4]) ;
    title('B. Inflation, \% Annual','Interpreter','latex') ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;    

print -depsc ./output/cf_fg_US_both_chipchiw_v4.eps
close  

tmp=100*(log(mean(output_chk'))-log(squeeze(mean(output_fg_US(:,ind_plot),2)')));
mean(tmp(find(data_dates==2009):find(data_dates==2015.5)))

tmp= 4*x_T(SET.variable.infl_obs,:)-4*squeeze(mean(y_fg_US_mat(4,:,ind_plot),3));
mean(tmp(find(data_dates==2009):find(data_dates==2015.5)))

%% Durations plots

figure ;

set(gcf, 'DefaultAxesLineWidth', 1);
set(gcf, 'DefaultLineLineWidth', 1.5);
set(gcf, 'DefaultAxesTickLabelInterpreter','latex'); 
set(gcf, 'DefaultLegendInterpreter','latex');

subplot(2,2,1) ; hold on ;
    plot(2009:0.25:2015.5,mean(params_T_f_),'-k','LineWidth',2)  
        plot(2009:0.25:2015.5,prctile(params_T_f_,[5 95]),':k','LineWidth',1)  
    line([2008 2016],[0 0],'Color','k','LineWidth',1) ;
    ylim([-2 20]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    xlim([2009 2016]) ;
    title('A. US, Actual Duration','Interpreter','latex') ;
    
subplot(2,2,2) ; hold on ;
	plot(2009:0.25:2015.5,mean(lb_T_vec(:,1:end)),'LineWidth',2,'Color',[0    0.4470    0.7410])
        plot(2009:0.25:2015.5,prctile(lb_T_vec,[5 95]),':','LineWidth',1,'Color',[0    0.4470    0.7410])  
    line([2008 2016],[0 0],'Color','k','LineWidth',1) ;
    ylim([-2 20]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    xlim([2009 2016]) ;
    title('B. US, Lower Bound Duration','Interpreter','latex') ;
    
subplot(2,2,3) ; hold on ;
    plot(2009.25:0.25:2014.75,mean(params_T_d_),'-k','LineWidth',2)
    plot(2009.25:0.25:2014.75,prctile(params_T_d_,[5 95]),':k','LineWidth',1)  
    line([2008 2016],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    title('C. Canada, Actual Duration','Interpreter','latex') ;
    xlim([2009 2016]) ;
    ylim([-2 20]) ;

subplot(2,2,4) ; hold on ;
    plot(2009.25:0.25:2014.75,mean(lb_T_d_vec(:,:)),'LineWidth',2,'Color',[0    0.4470    0.7410])  
    plot(2009.25:0.25:2014.75,prctile(lb_T_d_vec,[5 95]),':','LineWidth',1,'Color',[0    0.4470    0.7410])  
    line([2008 2016],[0 0],'Color','k','LineWidth',1) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
    title('D. Canada, Lower Bound Duration','Interpreter','latex') ;
    xlim([2009 2016]) ;
    ylim([-2 20]) ;

print -depsc ./output/durations_chipchiw_v4.eps
close     

%% Construct FG shocks
 
clear X2 X3 X4
 
X2 = [zeros(max_id,1), params_T_d_(ind,:)] ;
X3 = (params_T_f_(ind,:)-lb_T_vec) ;
X4 = (params_T_d_(ind,:)-lb_T_d_vec) ;

for jj=1:max_id
    fgshockUS(jj,1) = X3(jj,1) ;
    fgshocksCan(jj,1) = X4(jj,1) ;
    for ii=2:27
        if params_T_f_(ind(jj),ii-1)>0
            %disp('there')
            fgshockUS(jj,ii) = X3(jj,ii) - X3(jj,ii-1) ;
        end
        if params_T_f_(ind(jj),ii-1)==0
            %disp('hi')
            fgshockUS(jj,ii) = X3(jj,ii) - X3(jj,ii-1) + 1 ;
        end
        
        if ii<24
        if params_T_d_(ind(jj),ii-1)>0
            fgshocksCan(jj,ii) = X4(jj,ii) - X4(jj,ii-1) ;
        end
        if params_T_d_(ind(jj),ii-1)==0
            fgshocksCan(jj,ii) = X4(jj,ii) - X4(jj,ii-1) + 1 ;
        end       
        end
    end
end

for ii=1:23
tmp=corrcoef(fgshockUS(:,ii+1),fgshocksCan(:,ii));
tmp2(ii) = tmp(1,2) ;
end

CanFGDates=2009.25:0.25:2014.75;
USFGDates = 2009:.25:2015.5;

tmp3=corrcoef(vec(fgshockUS(:,find(USFGDates==2009.25):find(USFGDates==2014.75))),...
    vec(fgshocksCan(:,find(CanFGDates==2009.25):find(CanFGDates==2014.75))));
tmp3=tmp3(1,2) ;

tmp3=corrcoef(mean(fgshockUS(:,find(USFGDates==2009.25):find(USFGDates==2014.75))),...
    mean(fgshocksCan(:,find(CanFGDates==2009.25):find(CanFGDates==2014.75))));
tmp3=tmp3(1,2) ;

tmp4=corrcoef(mean(X3(:,find(USFGDates==2009.25):find(USFGDates==2014.75))),...
    mean(X4(:,find(CanFGDates==2009.25):find(CanFGDates==2014.75))));
tmp4=tmp4(1,2) ;

fprintf('\n')
fprintf('\n')
fprintf('Mode of FG duration:  US     = %9.2f \n', mode(X3(:)));
fprintf('Mode of FG duration:  Canada = %9.2f \n', mode(X4(:)));
fprintf('Mean of FG duration:  US     = %9.2f \n', mean(X3(:)));
fprintf('Mean of FG duration:  Canada = %9.2f \n', mean(X4(:)));
fprintf('Correlation of FG Shocks     = %9.2f \n', tmp3);
fprintf('Correlation of FG durations  = %9.2f \n', tmp4);


figure ;
set(gcf,'PaperPosition',[0 0 11 4]) ;
set(gcf,'DefaultLineLineWidth', 3) ;
subplot(1,2,1) ; hold on ;
    %[x,y]=hist(fgshockUS(:,10),13);%max(fgshockUS(:,10))-min(fgshockUS(:,10))) ;
    %bar(y,x/length(fgshockUS),1.0,'FaceColor',[0 0.4470 0.7410],'EdgeColor','k') ;
    ksdensity(fgshockUS(:,10),'NumPoints',50);
    line([0 0],[0 0.25],'Color','k','LineWidth',1,'LineStyle','--') ;
    xlim([-10 12]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',14,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
    xlabel('Quarters') ;
    title('A. 2011Q2') ;

subplot(1,2,2) ; hold on ;
    %[x,y]=hist(fgshockUS(:,18),12);%max(fgshockUS(:,18))-min(fgshockUS(:,18))) ;
   % bar(y,x/length(fgshockUS),1.0,'FaceColor',[0 0.4470 0.7410],'EdgeColor','k') ;
	%xlim([-20 20]) ;
    ksdensity(fgshockUS(:,18),'NumPoints',50) ;
    line([0 0],[0 0.25],'Color','k','LineWidth',1,'LineStyle','--') ;
    xlim([-14 10]) ;
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',14,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
    xlabel('Quarters') ;
    title('B. 2013Q2') ;  
    
set(findall(gcf, 'Type', 'Text'),'FontWeight', 'Normal')
print -depsc ./output/fig_fgshock_US2011q2_2013q2_chipchiw.eps
close


