function out = mcmc(x0, zlb_t, c, SET, Zobs, runs, cont_xx, cont_TT_d, cont_TT_f)
%out = mcmc(SET, data, zlb_t, hess, min_param_ll, runs)
% 
% Runs the Markov Chain Monte Carlo
%
% Outputs the chains for parameters and durations

if nargin < 6 ; runs = 1 ; end

%% Parameters

nparam = length(x0) ;

N      = SET.EST.N ; 
df     = SET.EST.df ; 

upperT = SET.EST.upperT ; % maximum number of updated elements in Tbs
Tbstar = SET.EST.Tbstar ; % Maximum expected duration at the ZLB is Tbstar;

pop = 1:nparam ;

kk_Tbs      = 0 ;
kk_xx       = 0 ; 

zlb_t_d = zlb_t(1,:) ;
zlb_t_f = zlb_t(2,:) ;

reject_rate = 0;

%% Posterior

load prior_Tbs_US.mat

PP_f = PP;

% Starts duration in the US at the mean of the prior
inidur_f = 1:size(PP_f,1); 
inidur_f = round(inidur_f*PP_f);
inidur_f = [inidur_f, 3, 2 ,1] ;

[rows_PP_f, cols_PP_f] = size(PP_f) ;

%if Tbstar < 23 ; disp('Increase Tbstar!') ; end

PP_f = 0.5*PP_f;
PP_f = [PP_f ; zeros(Tbstar,cols_PP_f)] ;
PP_f(Tbstar+1:end,:) = [];
PP_f = PP_f + 0.5/Tbstar;
PP_f = [PP_f, PP_f(:,end-2:end)] ; 

load prior_Tbs_canada.mat

PP_d = PP_can ;
[rows_PP_d, cols_PP_d] = size(PP_can) ;

PP_d = [PP_d ; zeros(Tbstar,cols_PP_d)] ;
PP_d(Tbstar+1:end,:) = []; 
PP_d(:,6:end) = 1/Tbstar ;
PP_d = 0.5*PP_d + 0.5/Tbstar;

posterior = @(x, T_d, T_f) ...
    - loglike(x, T_d, T_f, SET, Zobs) - prior(x) ...
    - log(prod(diag(PP_f(T_f(zlb_t_f>0),1:sum(zlb_t_f))))) ...
    - log(prod(diag(PP_d(T_d(zlb_t_d>0),1:sum(zlb_t_d))))) ;

%% Settings for durations and parameters

D_d = sum(zlb_t_d) ; 
D_f = sum(zlb_t_f) ; 

popTbs_d = 1:D_d ;
popTbs_f = 1:D_f ;

Tbs_prop_d = zeros(N,D_d); % Proposals
Tbs_prop_f = zeros(N,D_f); % Proposals

Tbs_prop_d(1,:) = 3 + Tbs_prop_d(1,:) ;
Tbs_prop_f(1,:) = inidur_f ;

Tbs_out_d = zeros(N,D_d);
Tbs_out_f = zeros(N,D_f);
T_d       = zeros(1,SET.ss) ;
T_f       = zeros(1,SET.ss) ;

Prop_xx      = zeros(N,length(x0));
Prop_xx(1,:) = x0';
xx_out       = zeros(N,length(x0));

Tbs_out_d(1,:) = Tbs_prop_d(1,:); % Initialize
Tbs_out_f(1,:) = Tbs_prop_f(1,:);
xx_out(1,:)    = x0';

%% If continuing a chain

if nargin>6 
    Prop_xx(1,:)    = cont_xx(end,:,runs) ;
    xx_out(1,:)     = cont_xx(end,:,runs) ;
    Tbs_prop_d(1,:) = cont_TT_d(end,:,runs) ;
    Tbs_out_d(1,:)  = cont_TT_d(end,:,runs) ;
    Tbs_prop_f(1,:) = cont_TT_f(end,:,runs) ;
    Tbs_out_f(1,:)  = cont_TT_f(end,:,runs) ;
end

%% DUAL SINGLE-MOVE SAMPLER
for j=2:N
%% T-block    

% Durations are same length as data. So zlb_t>0 finds the periods in
% T_f where the zero lower bound should bind
T_d(zlb_t_d>0) = Tbs_out_d(j-1,:) ;
T_f(zlb_t_f>0) = Tbs_out_f(j-1,:) ;

denominator = posterior(xx_out(j-1,:), T_d, T_f) ;

Tbs_prop_d(j,:) = Tbs_out_d(j-1,:);
Tbs_prop_f(j,:) = Tbs_out_f(j-1,:);

blocksizeTbs_d = ceil(upperT*rand);
blockTbs_d     = randsample(popTbs_d,blocksizeTbs_d);
Tbs_prop_d(j,blockTbs_d) = ceil((Tbstar)*rand(1,blocksizeTbs_d));

blocksizeTbs_f = ceil(upperT*rand);
blockTbs_f     = randsample(popTbs_f,blocksizeTbs_f);
Tbs_prop_f(j,blockTbs_f) = ceil((Tbstar)*rand(1,blocksizeTbs_f));

T_d(zlb_t_d>0) = Tbs_prop_d(j,:) ;
T_f(zlb_t_f>0) = Tbs_prop_f(j,:) ;

numerator = posterior(xx_out(j-1,:), T_d, T_f) ;

ratio = exp(-numerator + denominator) ;  
if isinf(ratio) ; ratio = 0 ; end
alpha = min(ratio,1) ;

if isnan(ratio) ; alpha = 0 ; end

if rand <= alpha
    Tbs_out_d(j,:) = Tbs_prop_d(j,:) ;
    Tbs_out_f(j,:) = Tbs_prop_f(j,:) ;
    kk_Tbs = kk_Tbs + 1;
    denominator2 = numerator ;
else
    Tbs_out_d(j,:) = Tbs_out_d(j-1,:);
    Tbs_out_f(j,:) = Tbs_out_f(j-1,:);
    denominator2 = denominator ;
end

acc_Tbs = 100*kk_Tbs/j;

%% xx-block

xx_prop     = xx_out(j-1,:);
nparam_max  = min(20,length(xx_prop));
blocksize   = ceil(nparam_max*rand);
block       = randsample(pop,blocksize);

xx_prop(1,block) = xx_out(j-1,block) + ...
    (c(block,block)*trnd(df,blocksize,1))'; % multivariate-t with df degrees of freedom

Prop_xx(j,:) = xx_prop(1,:);

T_d(zlb_t_d>0) = Tbs_out_d(j,:) ;
T_f(zlb_t_f>0) = Tbs_out_f(j,:) ;

% NOTE: imposes 0 < rhos < 1
if  min(xx_prop([1:4, 8:12,  21, 23:27, 31:38, 49:52])) <= 0 || max(xx_prop([1:4, 8:12,  21, 23:27, 31:38, 49:52])) >= 1  
    alpha2 = 0;
    ratio2 = 0;
    reject_rate = reject_rate + 1;
else
    numerator2   = posterior(xx_prop, T_d, T_f) ;
    ratio2 = exp(-numerator2 + denominator2);
    alpha2 = min(ratio2,1);
end

if isnan(ratio2) ; alpha2 = 0; end

if rand <= alpha2
   xx_out(j,:)=xx_prop;
   kk_xx = kk_xx + 1;
else
    xx_out(j,:) = xx_out(j-1,:);
end

acc_xx= 100*kk_xx/j;

if mod(j/100,1)==0 
    fprintf(sprintf('Chain %1.0f .... %2.3f per cent done\n      ***** accept rates: x, %2.3f and T, %2.3f\n',runs,100*j/N,acc_xx,acc_Tbs)) ;
end

if mod(j/100,1)==0 
    disp('mean of xx')
    mean(xx_out(1:j,:))
    disp('mode of T_d')
    mode(Tbs_out_d(1:j,:))
    disp('mode of T_f')
    mode(Tbs_out_f(1:j,:))
end

if mod(j/(N/SET.EST.numbersave),1)==0
    eval(['save ./output/chain_', int2str(runs), ';']) ;
end

end
%% Outputs

out.xx_out    = xx_out ;
out.Tbs_out_d = Tbs_out_d ;
out.Tbs_out_f = Tbs_out_f ;
out.acc_rates = [acc_xx, acc_Tbs, N, reject_rate] ;

