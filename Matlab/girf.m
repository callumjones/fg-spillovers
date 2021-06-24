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

%% Set estimated parameters and durations

load mhall_05-Mar-2020_estimatechis.mat

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

%% Select draws

date_for_irf = 2011.5 ;
FG_ext = 2 ;

tind = 2009:.25:2015.5;
pick_q = find(tind==date_for_irf);
pick_dur = 4 ; % Pick durations of length pick_dur

id = find(params_T_f_(:,pick_q)==pick_dur); 

if isempty(id)
    disp('No draw satisfies criteria')
end

s = RandStream('mlfg6331_64','Seed',17);

ndraw = 10000 ;
id    = id(randsample(s,1:length(id),ndraw)) ;
draws = id ;

maxdraw = length(draws) ; %change later to id_low also

%% Get shocks and run GIRFs

ss      = SET.ss ; 
n       = SET.variable.n_ ; 
nparam  = SET.nparam ;
nobs    = SET.nobs ;
no_T    = 0 ;

SC.t_f = zeros(1,SET.ss) ;
SC.t_d = zeros(1,SET.ss) ;    

SET.horizon = SET.ss ;

date_for_irf_idx = find(data_dates==date_for_irf) ;

for ii=1:maxdraw
        
%% 1. First extract the shocks for given posterior draw

SET.M_.params(SET.EST.params_to_estimate_idx) = params_x_(draws(ii),:) ;
SC.t_f(zlb_t_f>0) = params_T_f_(draws(ii),:) ;
SC.t_d(zlb_t_d>0) = params_T_d_(draws(ii),:) ;
dursest_f = params_T_f_(draws(ii),:) ;
dursest_d = params_T_d_(draws(ii),:) ;

dyn_in.M_       = SET.M_ ;
dyn_in.oo_      = SET.oo_ ;
dyn_in.options_ = SET.options_ ;
dyn_in.solve    = 1 ;

mats = resolve(SET,dyn_in) ; % Construct structural matrices: A0, A1, etc

SET.mats.Q = mats.mats.Q ;
SET.mats.G = mats.mats.G ;

SET.mat_init     = mats.mat_init ;
SET.mat_i_f_zlb  = mats.mat_i_f_zlb ;
SET.mat_i_d_zlb  = mats.mat_i_d_zlb ;
SET.mat_i_df_zlb = mats.mat_i_df_zlb ;
SET.mat_fin      = mats.mat_fin ;

mats = tv_mats(SET,SC) ; % Constructs Qt and Gt

Qt = mats.Qt ;
Gt = mats.Gt ;
    
Ct = Qt(1:n_-1,end,:) ;
Pt = Qt(1:n_-1,1:n_-1,:) ;
Dt = Gt(1:n_-1,:,:) ;

H = zeros(nobs,n_-1) ; 

for i=1:nobs
    H(i,SET.EST.obs(i)) = 1 ;
end

Omega     = eye(SET.variable.l_) ;
SET.Omega = Omega ;
    
[x_T, n_T] = ks_(SET, Ct, Pt, Dt, H, data, SC) ;

x_T = [x_T ; ones(1,ss) ] ;
sh = n_T ;
sh(:,(date_for_irf_idx+1):end) = 0 ;

%% Conventional IRF in 2006Q1

nosh = n_T ;
nosh(:,find(data_dates==2006):end) = 0 ;
sh_conv = nosh ; 
% sh_conv(SET.variable.shock.eps_r_f,find(data_dates==2006)) = -2.5 ; % gives similar size response to 2Q FG
sh_conv(SET.variable.shock.eps_r_f,find(data_dates==2006)) = -.65 ; 

y_sh_2006Q1{ii}   = x_T(:,1) ;
y_nosh_2006Q1{ii} = x_T(:,1) ;

for t=2:SET.ss
    y_sh_2006Q1{ii}(:,t) = Qt(:,:,1)*y_sh_2006Q1{ii}(:,t-1) + Gt(:,:,1)*sh_conv(:,t) ;
    y_nosh_2006Q1{ii}(:,t) = Qt(:,:,1)*y_nosh_2006Q1{ii}(:,t-1) + Gt(:,:,1)*nosh(:,t) ;
end

%% GIRF

% 1. Compute response without forward guidance shock

 % 
    SC.t_f = zeros(1,SET.ss) ;
    SC.t_d = zeros(1,SET.ss) ;

    SC.t_f(zlb_t_f>0) = dursest_f ;
    SC.t_f(date_for_irf_idx) = SC.t_f(date_for_irf_idx);
    for tt=(date_for_irf_idx+1):SET.ss, SC.t_f(tt) = SC.t_f(tt-1)-1 ; end
    SC.t_f(SC.t_f<0) = 0 ;

    SC.t_d(zlb_t_d>0) = dursest_d ;
    for tt=(date_for_irf_idx+1):SET.ss, SC.t_d(tt) = SC.t_d(tt-1)-1 ; end
    SC.t_d(SC.t_d<0) = 0 ;
    
    mats = tv_mats(SET,SC) ;

    Qt = mats.Qt ;
    Gt = mats.Gt ; 

    Qt(:,:,SET.ss+1:(SET.ss+SET.maxchk)) = ...
       repmat(Qt(:,:,end),[1 1 SET.maxchk]) ;
    Gt(:,:,SET.ss+1:(SET.ss+SET.maxchk)) = ...
       repmat(Gt(:,:,end),[1 1 SET.maxchk]) ;

     y_fg{ii} = x_T(:,1) ;

     for t=2:SET.ss
         y_fg{ii}(:,t) = Qt(:,:,t)*y_fg{ii}(:,t-1) +  Gt(:,:,t)*sh(:,t) ;
     end

% 2. Next compute the response with forward guidance shock

    % More US FG by 2 quarters
    SC.t_f = zeros(1,SET.ss) ;
    SC.t_d = zeros(1,SET.ss) ;

    SC.t_f(zlb_t_f>0) = dursest_f ;
    SC.t_f(date_for_irf_idx) = SC.t_f(date_for_irf_idx) + FG_ext;
    for tt=(date_for_irf_idx+1):SET.ss, SC.t_f(tt) = SC.t_f(tt-1)-1 ; end
    SC.t_f(SC.t_f<0) = 0 ;

    SC.t_d(zlb_t_d>0) = dursest_d ;
    for tt=(date_for_irf_idx+1):SET.ss, SC.t_d(tt) = SC.t_d(tt-1)-1 ; end
    SC.t_d(SC.t_d<0) = 0 ;
    
    mats = tv_mats(SET,SC) ;

    Qt = mats.Qt ;
    Gt = mats.Gt ; 

    Qt(:,:,SET.ss+1:(SET.ss+SET.maxchk)) = ...
        repmat(Qt(:,:,end),[1 1 SET.maxchk]) ;
    Gt(:,:,SET.ss+1:(SET.ss+SET.maxchk)) = ...
        repmat(Gt(:,:,end),[1 1 SET.maxchk]) ;
    
    y_fg_ext{ii} = x_T(:,1) ;

    for t=2:SET.ss
        y_fg_ext{ii}(:,t) = Qt(:,:,t)*y_fg_ext{ii}(:,t-1) +  Gt(:,:,t)*sh(:,t) ;
    end
    
fprintf('\n')
fprintf('Draw %i, US forward guidance extension of %i quarters from duration %i in %4.1f \n', ii, FG_ext, pick_dur, date_for_irf)

end

%% Plot

y_sh   = zeros(n,ss,maxdraw);
y_nosh = zeros(n,ss,maxdraw);
y_     = zeros(n,ss,maxdraw);
y_ext  = zeros(n,ss,maxdraw);

for i_ = 1:maxdraw 
    y_sh(:,:,i_)   = y_sh_2006Q1{i_} ; 
    y_nosh(:,:,i_) = y_nosh_2006Q1{i_};
    y_(:,:,i_)     = y_fg{i_} ;
    y_ext(:,:,i_)  = y_fg_ext{i_} ;
end

y_con   = y_sh - y_nosh ; 
y_uncon = y_ext - y_ ; 

st  = find(data_dates==2006) ;
st2 = find(data_dates==date_for_irf) ;

lw1 = 0.25;
lw2 = 2;

up_to = 16;

%%

SwatheOpt = PlotSwatheOption;
SwatheOpt.marker = '*';
SwatheOpt.trans = 1;
SwatheOpt.xaxis = (1:up_to);
SwatheOpt.swatheonly = 1 ;

figure; 
set(gcf,'DefaultLineLineWidth', 2) ;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 4]) ;
subplot(2,3,1) ; hold on ;

    plot(1:up_to,mean(100*squeeze((y_con(SET.variable.y_f,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(100*squeeze(y_uncon(SET.variable.y_f,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 
    
    SwatheOpt.swathecol  = [150, 150, 150]./255; 

    PlotSwathe(100*squeeze(prctile(y_con(SET.variable.y_f,st:st+up_to-1,:),50,3)),...
        100*squeeze(prctile(y_con(SET.variable.y_f,st:st+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    SwatheOpt.swathecol  = [138, 178, 212]./255; 

    PlotSwathe(100*squeeze(prctile(y_uncon(SET.variable.y_f,st2:st2+up_to-1,:),50,3)),...
        100*squeeze(prctile(y_uncon(SET.variable.y_f,st2:st2+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    
    plot(1:up_to,mean(100*squeeze((y_con(SET.variable.y_f,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(100*squeeze(y_uncon(SET.variable.y_f,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 

    line([0 up_to],[0 0],'Color','k','LineWidth',1) ;
    xlim([0 up_to]) ;
    ylim([0 0.6]) ;
    title('A. US Output, \%','Interpreter','latex') ;
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')     
            
    set(gca,'xtick',[0:4:up_to])
    
    h=legend('Conventional','2Q Fwd Guidance','location','northeast');
    set(h,'Interpreter','latex','FontSize',6);

subplot(2,3,2) ; hold on ;

    plot(1:up_to,mean(1*squeeze((y_con(SET.variable.infl_f_obs,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(1*squeeze(y_uncon(SET.variable.infl_f_obs,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 
    
	SwatheOpt.swathecol  = [150, 150, 150]./255; 

    PlotSwathe(1*squeeze(prctile(y_con(SET.variable.infl_f_obs,st:st+up_to-1,:),50,3)),...
        1*squeeze(prctile(y_con(SET.variable.infl_f_obs,st:st+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    SwatheOpt.swathecol  = [138, 178, 212]./255; 

    PlotSwathe(1*squeeze(prctile(y_uncon(SET.variable.infl_f_obs,st2:st2+up_to-1,:),50,3)),...
        1*squeeze(prctile(y_uncon(SET.variable.infl_f_obs,st2:st2+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    
    plot(1:up_to,mean(1*squeeze((y_con(SET.variable.infl_f_obs,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(1*squeeze(y_uncon(SET.variable.infl_f_obs,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 

    line([0 up_to],[0 0],'Color','k','LineWidth',1) ;
    xlim([0 up_to]) ;
    title('B. US Inflation, \%','Interpreter','latex') ;

    ylim([0 0.08]) ;
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')    
    
    set(gca,'xtick',[0:4:up_to])
    
subplot(2,3,3) ; hold on ;

    plot(1:up_to,mean(1*squeeze((y_con(SET.variable.r_f_1_obs,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(1*squeeze(y_uncon(SET.variable.r_f_1_obs,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 
    
	SwatheOpt.swathecol  = [150, 150, 150]./255; 

    PlotSwathe(1*squeeze(prctile(y_con(SET.variable.r_f_1_obs,st:st+up_to-1,:),50,3)),...
        1*squeeze(prctile(y_con(SET.variable.r_f_1_obs,st:st+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    SwatheOpt.swathecol  = [138, 178, 212]./255; 

    PlotSwathe(1*squeeze(prctile(y_uncon(SET.variable.r_f_1_obs,st2:st2+up_to-1,:),50,3)),...
        1*squeeze(prctile(y_uncon(SET.variable.r_f_1_obs,st2:st2+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    
    plot(1:up_to,mean(1*squeeze((y_con(SET.variable.r_f_1_obs,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(1*squeeze(y_uncon(SET.variable.r_f_1_obs,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 
 
    line([0 up_to],[0 0],'Color','k','LineWidth',1) ;
    xlim([0 up_to]) ;
    title('C. Fed Funds Rate, Annual \%','Interpreter','latex') ;

    ylim([-0.4 0]) ;
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')    
    
    set(gca,'xtick',[0:4:up_to])

subplot(2,3,4) ; hold on ;

    plot(1:up_to,mean(100*squeeze((y_con(SET.variable.y,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(100*squeeze(y_uncon(SET.variable.y,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 
 
    SwatheOpt.swathecol  = [150, 150, 150]./255; 

    PlotSwathe(100*squeeze(prctile(y_con(SET.variable.y,st:st+up_to-1,:),50,3)),...
        100*squeeze(prctile(y_con(SET.variable.y,st:st+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    SwatheOpt.swathecol  = [138, 178, 212]./255; 

    PlotSwathe(100*squeeze(prctile(y_uncon(SET.variable.y,st2:st2+up_to-1,:),50,3)),...
        100*squeeze(prctile(y_uncon(SET.variable.y,st2:st2+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    
    plot(1:up_to,mean(100*squeeze((y_con(SET.variable.y,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(100*squeeze(y_uncon(SET.variable.y,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 

    line([0 up_to],[0 0],'Color','k','LineWidth',1) ;
    xlim([0 up_to]) ;
    title('D. Canada Output, \%','Interpreter','latex') ;
    ylim([-0.6 0]) ;
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')    
    
    set(gca,'xtick',[0:4:up_to])

subplot(2,3,5) ; hold on ;

    plot(1:up_to,mean(1*squeeze((y_con(SET.variable.infl_obs,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(1*squeeze(y_uncon(SET.variable.infl_obs,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 
    
    SwatheOpt.swathecol  = [150, 150, 150]./255; 

    PlotSwathe(1*squeeze(prctile(y_con(SET.variable.infl_obs,st:st+up_to-1,:),50,3)),...
        1*squeeze(prctile(y_con(SET.variable.infl_obs,st:st+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    SwatheOpt.swathecol  = [138, 178, 212]./255; 

    PlotSwathe(1*squeeze(prctile(y_uncon(SET.variable.infl_obs,st2:st2+up_to-1,:),50,3)),...
        1*squeeze(prctile(y_uncon(SET.variable.infl_obs,st2:st2+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    
    plot(1:up_to,mean(1*squeeze((y_con(SET.variable.infl_obs,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(1*squeeze(y_uncon(SET.variable.infl_obs,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 

    line([0 up_to],[0 0],'Color','k','LineWidth',1) ;
    xlim([0 up_to]) ;
    title('E. Canada Inflation, \%','Interpreter','latex') ;
    %ylabel('%') ;
    ylim([-0.15 0.05000001]) ;
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')   
    
    set(gca,'xtick',[0:4:up_to])

subplot(2,3,6) ; hold on ;

    plot(1:up_to,mean(100*squeeze((y_con(SET.variable.q,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(100*squeeze(y_uncon(SET.variable.q,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 
   
    SwatheOpt.swathecol  = [150, 150, 150]./255; 

    PlotSwathe(100*squeeze(prctile(y_con(SET.variable.q,st:st+up_to-1,:),50,3)),...
        100*squeeze(prctile(y_con(SET.variable.q,st:st+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    SwatheOpt.swathecol  = [138, 178, 212]./255; 

    PlotSwathe(100*squeeze(prctile(y_uncon(SET.variable.q,st2:st2+up_to-1,:),50,3)),...
        100*squeeze(prctile(y_uncon(SET.variable.q,st2:st2+up_to-1,:),[5 95],3)),SwatheOpt); hold on;
    
    plot(1:up_to,mean(100*squeeze((y_con(SET.variable.q,st:st+up_to-1,:)))'),'-k','LineWidth',lw2)    
    plot(1:up_to,mean(100*squeeze(y_uncon(SET.variable.q,st2:st2+up_to-1,:))'),'Color',[0 0.4470 0.7410],'LineWidth',lw2) 

    line([0 up_to],[0 0],'Color','k','LineWidth',1) ;
    xlim([0 up_to]) ;
    title('F. Real Exchange Rate, \%','Interpreter','latex') ;
    ylim([-1.2 0]) ;
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')    
    
    set(gca,'xtick',[0:4:up_to])
   
set(findall(gcf, 'Type', 'Text'),'FontWeight', 'Normal')
print -depsc ./output/fg_girf_2Q_4_2011q3.eps
close
