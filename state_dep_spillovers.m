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

date_for_irf = 2013.5 ; % 2013Q3
FG_ext = 2 ;

tind = 2009:.25:2014.75;
pick_q = find(tind==date_for_irf);
pick_dur = 8 ; % Pick durations of length pick_dur
highxif = 0 ;
lowxif  = 0 ;
meanxif = 1 ;

if highxif
    id = find(and(params_T_f_(:,pick_q)==pick_dur, ...
        params_x_(:,find(strcmp(SET.EST.params_to_estimate,'sig_xi_f')))>0.4)); % hixif 
elseif lowxif
    id = find(and(params_T_f_(:,pick_q)==pick_dur, ...
        params_x_(:,find(strcmp(SET.EST.params_to_estimate,'sig_xi_f')))<0.18)); % lowxif 
elseif meanxif 
    id = find(params_T_f_(:,pick_q)==pick_dur); % meanxif 
end    

if isempty(id)
    disp('No draw satisfies criteria')
    return
end

ndraw = 100 ;
id    = id(randsample(1:length(id),ndraw)) ;
draws = id ;

maxdraw = length(draws) ; 

%% Get shocks and run GIRFs

ss      = SET.ss ; 
n       = SET.variable.n_ ; 
nparam  = SET.nparam ;
nobs    = SET.nobs ;
no_T    = 0 ;

SC.t_f = zeros(1,SET.ss) ;
SC.t_d = zeros(1,SET.ss) ;    

SET.horizon = SET.ss ;

for ii=1:maxdraw
        
% 1. First extract the shocks for given posterior draws

SET.M_.params(SET.EST.params_to_estimate_idx) = params_x_(draws(ii),:) ; 
SC.t_f(zlb_t_f>0) = params_T_f_(draws(ii),:) ; 
SC.t_d(zlb_t_d>0) = params_T_d_(draws(ii),:) ; 

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

y_zlb     = x_T(:,1) ;
z_f_T     = zeros(1,SET.ss) ;
z_d_T     = zeros(1,SET.ss) ;

Qt_nozlb = repmat(Qt(:,:,1),[1 1 SET.ss+SET.maxchk]) ;
Gt_nozlb = repmat(Gt(:,:,1),[1 1 SET.ss+SET.maxchk]) ;

for t=2:SET.ss
    [y_zlb(:,t), ~, ~, z_f_T(t), z_d_T(t)] = ...
        zlb(SET, y_zlb(:,t-1), n_T(:,t), Qt_nozlb, Gt_nozlb, t) ;
end


% 2. Next compute the response when there is no forward guidance in Canada

date_for_irf_idx = find(data_dates==date_for_irf) ;

SC.t_f = zeros(1,SET.ss) ;
SC.t_d = zeros(1,SET.ss) ;

SC.t_f(zlb_t_f>0) = params_T_f_(draws(ii),:) ; 
for tt=(date_for_irf_idx+1):SET.ss ; SC.t_f(tt) = SC.t_f(tt-1)-1 ; end
SC.t_f(SC.t_f<0) = 0 ;

SC.t_d(zlb_t_d>0) = 0*params_T_d_(draws(ii),:) ;
SC.t_d(date_for_irf_idx:end) = 0 ;

y_no_fg{ii} = x_T(:,1) ;

sh = n_T ;
sh(:,(date_for_irf_idx+1):end) = 0 ;

mats = tv_mats(SET,SC) ;

Qt = mats.Qt ;
Gt = mats.Gt ; 

Qt(:,:,SET.ss+1:(SET.ss+SET.maxchk)) = ...
    repmat(Qt(:,:,end),[1 1 SET.maxchk]) ;
Gt(:,:,SET.ss+1:(SET.ss+SET.maxchk)) = ...
    repmat(Gt(:,:,end),[1 1 SET.maxchk]) ;

for t=2:SET.ss
    y_no_fg{ii}(:,t) = Qt(:,:,t)*y_no_fg{ii}(:,t-1) +  Gt(:,:,t)*sh(:,t) ;
end


% 3. Next compute the response when there is no FG in Canada and extention 
% in the US

SC.t_f = zeros(1,SET.ss) ;
SC.t_d = zeros(1,SET.ss) ;

SC.t_f(zlb_t_f>0) = params_T_f_(draws(ii),:) ; 

SC.t_f(date_for_irf_idx) = SC.t_f(date_for_irf_idx)+1*FG_ext;

for tt=(date_for_irf_idx+1):SET.ss ; SC.t_f(tt) = SC.t_f(tt-1)-1 ; end
SC.t_f(SC.t_f<0) = 0 ;

SC.t_d(zlb_t_d>0) = 0*params_T_d_(draws(ii),:) ;
SC.t_d(date_for_irf_idx:end) = 0 ;

y_fg_ext{ii} = x_T(:,1) ;

mats = tv_mats(SET,SC) ;

Qt = mats.Qt ;
Gt = mats.Gt ; 

Qt(:,:,SET.ss+1:(SET.ss+SET.maxchk)) = ...
    repmat(Qt(:,:,end),[1 1 SET.maxchk]) ;
Gt(:,:,SET.ss+1:(SET.ss+SET.maxchk)) = ...
    repmat(Gt(:,:,end),[1 1 SET.maxchk]) ;

for t=2:SET.ss
    y_fg_ext{ii}(:,t) = Qt(:,:,t)*y_fg_ext{ii}(:,t-1) +  Gt(:,:,t)*sh(:,t) ;
end

fprintf('\n')
fprintf('Draw %i, US forward guidance extension of %i quarters in %4.1f \n', ii, FG_ext, date_for_irf)

end

%% Generate GIRFs and series to plot

toplot = [] ;

for ii=1:length(y_fg_ext)
    toplot(:,:,ii)  = y_fg_ext{ii} - y_no_fg{ii} ; % GIRF
    xi_f_draw(ii,:) = y_fg_ext{ii}(SET.variable.xi_f,:) ;
    xi_draw(ii,:)   = y_fg_ext{ii}(SET.variable.xi,:) ;
    rp_draw(ii,:)   = y_fg_ext{ii}(SET.variable.rp,:) ;
    a_f_draw(ii,:)  = y_fg_ext{ii}(SET.variable.a_f,:) ;
    a_draw(ii,:)    = y_fg_ext{ii}(SET.variable.a,:) ;
end

[a,date_for_irf_idx_2] = sort(xi_f_draw) ;

toplot = toplot(:,:,date_for_irf_idx_2(:,1)) ;

if highxif
    
eval(['save ./output/girfs_state_dep_',int2str(date_for_irf),'_',...
    int2str(pick_dur),'_',int2str(FG_ext),'_hixif toplot xi_f_draw xi_draw data_dates']) ;

elseif lowxif
    
eval(['save ./output/girfs_state_dep_',int2str(date_for_irf),'_',...
    int2str(pick_dur),'_',int2str(FG_ext),'_loxif toplot xi_f_draw xi_draw data_dates']) ;

elseif meanxif

eval(['save ./output/girfs_state_dep_',int2str(date_for_irf),'_',...
    int2str(pick_dur),'_',int2str(FG_ext),'_meanxif toplot xi_f_draw xi_draw data_dates']) ;
  
end

%%

load ./output/girfs_state_dep_2014_8_2_hixif.mat

tmp = 1:4:size(toplot,3);

figure ;
set(gcf,'DefaultLineLineWidth', 2) ;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 4]) ;
subplot(2,3,1) ; hold on ;

    load ./output/girfs_state_dep_2014_8_2_hixif.mat
    
    plot(data_dates,100*mean(squeeze(toplot(SET.variable.y_f,1:length(data_dates),:)),2)','-.','Color',[0 0.4470 0.7410],'LineWidth',2);
      
    load ./output/girfs_state_dep_2014_8_2_loxif.mat
    
    plot(data_dates,100*mean(squeeze(toplot(SET.variable.y_f,1:length(data_dates),:)),2)','-','Color',[0 0 0],'LineWidth',2);
        
    line([2010 2018],[0 0],'Color','k','LineWidth',1) ;
    title('A. US Output, \%','Interpreter','latex')
    xlim([2013 2018]);
    ylim([0 0.80000001])
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')     
    
    set(gca,'DefaultTextInterpreter','latex','TickLabelInterpreter','latex')   
    
subplot(2,3,2) ; hold on ;

    load ./output/girfs_state_dep_2014_8_2_hixif.mat

    plot(data_dates,1*mean(squeeze(toplot(SET.variable.infl_f_obs,1:length(data_dates),:)),2)','-.','Color',[0 0.4470 0.7410],'LineWidth',2);
      
    load ./output/girfs_state_dep_2014_8_2_loxif.mat
    
    plot(data_dates,1*mean(squeeze(toplot(SET.variable.infl_f_obs,1:length(data_dates),:)),2)','-','Color',[0 0 0],'LineWidth',2);
        
    line([2010 2018],[0 0],'Color','k','LineWidth',1) ;
    title('B. US Inflation, \%','Interpreter','latex')
    xlim([2013 2018]);%,ylim([0 .25])   
    
    h=legend('Large $\xi^\star$', 'Small $\xi^\star$','location','northeast') ;
    set(h,'FontSize',6.0,'Interpreter','latex');
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')     
    
    set(gca,'DefaultTextInterpreter','latex','TickLabelInterpreter','latex') 

subplot(2,3,3) ; hold on ;

    load ./output/girfs_state_dep_2014_8_2_hixif.mat

    plot(data_dates,1*mean(squeeze(toplot(SET.variable.r_f_1_obs,1:length(data_dates),:)),2)','-.','Color',[0 0.4470 0.7410],'LineWidth',2);
      
    load ./output/girfs_state_dep_2014_8_2_loxif.mat
    
    plot(data_dates,1*mean(squeeze(toplot(SET.variable.r_f_1_obs,1:length(data_dates),:)),2)','-','Color',[0 0 0],'LineWidth',2);
        
    line([2010 2018],[0 0],'Color','k','LineWidth',1) ;
    title('C. Fed Funds Rate, Annual \%','Interpreter','latex')
    xlim([2013 2018]);
    ylim([-0.6 0.1]) ;
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')     
    
    set(gca,'DefaultTextInterpreter','latex','TickLabelInterpreter','latex') 

subplot(2,3,4) ; hold on ;

    load ./output/girfs_state_dep_2014_8_2_hixif.mat

    plot(data_dates,100*mean(squeeze(toplot(SET.variable.y,1:length(data_dates),:)),2)','-.','Color',[0 0.4470 0.7410],'LineWidth',2);
    
    load ./output/girfs_state_dep_2014_8_2_loxif.mat

    plot(data_dates,100*mean(squeeze(toplot(SET.variable.y,1:length(data_dates),:)),2)','-','Color',[0 0 0],'LineWidth',2);
        
    line([2010 2018],[0 0],'Color','k','LineWidth',1) ;
        
    title('D. Canada Output, \%','Interpreter','latex')
    xlim([2013 2018]);
    ylim([-0.3 0.1])
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')     
    
    set(gca,'DefaultTextInterpreter','latex','TickLabelInterpreter','latex') 
    
subplot(2,3,5) ; hold on ;

    load ./output/girfs_state_dep_2014_8_2_hixif.mat

    plot(data_dates,1*mean(squeeze(toplot(SET.variable.infl_obs,1:length(data_dates),:)),2)','-.','Color',[0 0.4470 0.7410],'LineWidth',2);
    
    load ./output/girfs_state_dep_2014_8_2_loxif.mat

    plot(data_dates,1*mean(squeeze(toplot(SET.variable.infl_obs,1:length(data_dates),:)),2)','-','Color',[0 0 0],'LineWidth',2);
        
    line([2010 2018],[0 0],'Color','k','LineWidth',1) ;
    title('E. Canada Inflation, \%','Interpreter','latex')
    xlim([2013 2018]);
    ylim([-0.06 0.01])   
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')     
    
    set(gca,'DefaultTextInterpreter','latex','TickLabelInterpreter','latex') 
    
subplot(2,3,6) ; hold on ;

    load ./output/girfs_state_dep_2014_8_2_hixif.mat

    plot(data_dates,100*mean(squeeze(toplot(SET.variable.q,1:length(data_dates),:)),2)','-.','Color',[0 0.4470 0.7410],'LineWidth',2);
    
    load ./output/girfs_state_dep_2014_8_2_loxif.mat

    plot(data_dates,100*mean(squeeze(toplot(SET.variable.q,1:length(data_dates),:)),2)','-','Color',[0 0 0],'LineWidth',2); %[0.8500    0.3250    0.0980]
        
    line([2010 2018],[0 0],'Color','k','LineWidth',1) ;
    title('F. Real Exchange Rate, \%','Interpreter','latex')
    xlim([2013 2018]);
    ylim([-1.2 0.1])   
    
    ax = gca ;
    ax.YGrid='on';
    box('off') ;
	set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10,'LineWidth',1.5)
    set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
    
    set(gca,'DefaultTextInterpreter','latex','TickLabelInterpreter','latex')
    
set(findall(gcf, 'Type', 'Text'),'FontWeight', 'Normal')

print -depsc ./output/fig_state_dep_xif_chipchiw_v2.eps

close ;       
