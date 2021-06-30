function out = draw_from_prior(SET, x0)
%out = mcmc(SET, data, zlb_t, c, x0, runs, cont_xx, cont_TT_d, cont_TT_f)
% 
% Runs the Markov Chain Monte Carlo for the prior
%

runs=1;

%% Parameters

nparam = length(x0) ;

N      = 300000 ; 

pop = 1:nparam ;

kk_xx       = 0 ;
reject_rate = 0 ; 

%% Settings for durations and parameters

Prop_xx      = zeros(N,length(x0));
Prop_xx(1,:) = x0';
xx_out       = zeros(N,length(x0));

xx_out(1,:)  = x0';

load c

df = 1 ;

SET.estimate_tr = 1;

denominator2 = -prior(xx_out(1,:)) ;

%% DUAL SINGLE-MOVE SAMPLER
for j=2:N

%% xx-block

xx_prop     = xx_out(j-1,:);
blocksize   = 8;
block       = randsample(pop,8);

xx_prop(1,block) = xx_out(j-1,block) + ...
    (2*c(block,block)*trnd(df,blocksize,1))'; % multivariate-t with df degrees of freedom

Prop_xx(j,:) = xx_prop(1,:);

if sum(Prop_xx(j,:)<0)>0 ; xx_out(j,:)=xx_out(j-1,:) ; continue ; end

numerator2   = -prior(xx_prop) ;

ratio2 = exp(-numerator2 + denominator2);
alpha2 = min(ratio2,1);

if isnan(ratio2) ; alpha2 = 0 ; end

if rand < alpha2
   xx_out(j,:)=xx_prop;
   kk_xx = kk_xx + 1;
   denominator2=numerator2;
else
    xx_out(j,:) = xx_out(j-1,:);
end

acc_xx= 100*kk_xx/j;

if mod(j/10,1)==0 
    fprintf(sprintf('Chain %1.0f .... %2.3f per cent done\n      ***** accept rates: x, %2.3f\n',runs,100*j/N,acc_xx)) ;
end

end
%% Outputs

out   = xx_out(1000:end,:) ;

