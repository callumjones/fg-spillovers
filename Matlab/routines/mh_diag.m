function [sq_R, B, W, var_overall] = mh_diag(theta)

theta = permute(theta,[2 1 3]) ;

num_vars = length(theta(:,1,1)) ;
num_runs = length(theta(1,1,:)) ;
num_iters = length(theta(1,:,1)) ;

mean_across_chains  = zeros(num_vars,num_runs) ;
s                   = zeros(num_vars,num_runs) ;
  
mean_across_chains(1:num_vars,1:num_runs)   = mean(theta,2);
grand_mean                                  = mean(mean(theta,3),2) ;
s(1:num_vars,1:num_runs)                    = var(theta,0,2) ;

B = (num_iters / (num_runs - 1) * ...
    sum(((mean_across_chains - grand_mean*ones(1,num_runs)).^2),2)) ;
W = mean(s,2) ;

var_overall = ((num_iters - 1) / (num_iters)) * W + (1 / num_iters) * B ;
sq_R        = sqrt(var_overall ./ W ) ;