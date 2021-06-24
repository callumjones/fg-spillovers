load mhall_05-Mar-2020_estimatechis.mat

num_params = length(SET.EST.params_to_estimate) ;

figure(1) ;
for ii=1:16
subplot(4,4,ii); hold on ;
plot(params_x(:,ii,1)) ;
plot(params_x(:,ii,2)) ;
title(SET.EST.params_to_estimate{ii}) ;
end

figure(2) ;
for ii=17:32
subplot(4,4,ii-16); hold on ;
plot(params_x(:,ii,1)) ;
plot(params_x(:,ii,2)) ;
title(SET.EST.params_to_estimate{ii}) ;    
end

figure(3) ;
for ii=33:48
subplot(4,4,ii-32); hold on ;
plot(params_x(:,ii,1)) ;
plot(params_x(:,ii,2)) ;
title(SET.EST.params_to_estimate{ii}) ;    
end
    
% Durations f

figure(4) ; 
for ii = 1:16
    subplot(4,4,ii), hist(params_T_f_(:,ii),12), title(ii)
end

figure(5) ; 
for ii = 17:27
    subplot(4,4,ii-16), hist(params_T_f_(:,ii),12), title(ii)
end

%% Plot prior and posteriors

load mhall_05-Mar-2020_estimatechis.mat

num_params = length(SET.EST.params_to_estimate) ;

%priorx = draw_from_prior(SET, params_x(end,:,1)) ; save prior_sav priorx
load prior_sav % if done before

params_x_   = [] ;
params_T_d_ = [] ;
params_T_f_ = [] ;

burn = 50000 ;

for runs=1:maxproc
    params_x_ = [params_x_ ; params_x(burn:end,:,runs)] ; 
    params_T_d_ = [params_T_d_ ; params_T_d(burn:end,:,runs)] ;
    params_T_f_ = [params_T_f_ ; params_T_f(burn:end,:,runs)] ; 
end

figure ;
for ii=1:9
subplot(3,3,ii); hold on ;
ksdensity(priorx(:,ii)) ;
ksdensity(params_x_(:,ii)) ;
title(SET.EST.params_to_estimate{ii}) ;
end
print -depsc ./output/prior_post_plot1.eps
close

figure ;
for ii=10:18
subplot(3,3,ii-9); hold on ;
ksdensity(priorx(:,ii)) ;
ksdensity(params_x_(:,ii)) ;
title(SET.EST.params_to_estimate{ii}) ;
end
print -depsc ./output/prior_post_plot2.eps
close

figure ;
for ii=19:27
subplot(3,3,ii-18); hold on ;
ksdensity(priorx(:,ii)) ;
ksdensity(params_x_(:,ii)) ;
title(SET.EST.params_to_estimate{ii}) ;
end
print -depsc ./output/prior_post_plot3.eps
close

figure ;
for ii=28:36
subplot(3,3,ii-27); hold on ;
ksdensity(priorx(:,ii)) ;
ksdensity(params_x_(:,ii)) ;
title(SET.EST.params_to_estimate{ii}) ;
end
print -depsc ./output/prior_post_plot4.eps
close

figure ;
for ii=37:45
subplot(3,3,ii-36); hold on ;
ksdensity(priorx(:,ii)) ;
ksdensity(params_x_(:,ii)) ;
title(SET.EST.params_to_estimate{ii}) ;
end
print -depsc ./output/prior_post_plot5.eps
close

figure ;
for ii=46:52
subplot(3,3,ii-45); hold on ;
%ksdensity(priorx(:,ii)) ;
ksdensity(params_x_(:,ii)) ;
title(SET.EST.params_to_estimate{ii}) ;
end
print -depsc ./output/prior_post_plot6.eps
close
