function [fprior, p] = prior(x)

% Prior for x

p = zeros(1,length(x));

j=1 ;

% h_f
p(j) = log(betapdf(x(j),58.1, 24.9)); j=j+1 ;

% theta_p_f
p(j) = log(betapdf(x(j),224.25, 74.75)); j=j+1 ;

% theta_w_f
p(j) = log(betapdf(x(j),224.25, 74.75)); j=j+1 ;

% rho_r_f
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% phi_pi_f
p(j) = log(normpdf(x(j),2,0.25)); j=j+1 ;

% phi_g_f 
p(j) = log(gampdf(x(j),16, 0.03125)); j=j+1 ;

% phi_y_f 
p(j) = log(gampdf(x(j),16, 0.03125)); j=j+1 ;

% rho_xi_f 
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% rho_g_f
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% rho_xi_p_f
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% rho_xi_w_f
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% rho_tp_f_8
p(j) =  log(betapdf(x(j),2.625,2.625)); j=j+1 ;

% sig_z
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_r_f
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_xi_f
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_g_f
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_xi_p_f
p(j) = log(invgam_pdf(x(j),6,0.75)); j = j + 1 ;

% sig_xi_w_f
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_eps_r_f_8
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% r_f_8_const
p(j) =  log(normpdf(x(j),0.0,0.5)); j=j+1 ; % 8q

% %%% DOMESTIC PARAMETERS

% h
p(j) = log(betapdf(x(j),58.1, 24.9)); j=j+1 ;

% tau
p(j) = log(normpdf(x(j),1,0.5)) ; j=j+1; 

% theta_p
p(j) = log(betapdf(x(j),224.25, 74.75)); j=j+1 ;

% theta_w
p(j) = log(betapdf(x(j),224.25, 74.75)); j=j+1 ;

% theta_X
p(j) = log(betapdf(x(j),224.25, 74.75)); j=j+1 ;  

% theta_F
p(j) = log(betapdf(x(j),224.25, 74.75)); j=j+1 ; 

% rho_r 
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% phi_pi
p(j) = log(normpdf(x(j),2,0.25)); j=j+1 ;

% phi_g  
p(j) = log(gampdf(x(j),16, 0.03125)); j=j+1 ;

% phi_y 
p(j) = log(gampdf(x(j),16, 0.03125)); j=j+1 ;

% rho_rp
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% rho_xi
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% rho_g
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;
 
% rho_xi_H
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% rho_xi_w
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% rho_xi_X
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;
 
% rho_xi_F
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;
 
% rho_tp_8
p(j) = log(betapdf(x(j),2.625, 2.625)); j=j+1 ;

% sig_r
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_rp 
p(j) = log(invgam_pdf(x(j),6,0.75)); j = j + 1 ;

% sig_g 
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_xi
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_xi_H
p(j) = log(invgam_pdf(x(j),13.1,3.02)); j = j + 1 ;

% sig_xi_w
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;

% sig_xi_X
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;
 
% sig_xi_F
p(j) = log(invgam_pdf(x(j),6,0.75)); j = j + 1 ;
    
% sig_eps_r_8
p(j) = log(invgam_pdf(x(j),2.0051,0.50255)); j = j + 1 ;
 
% r_8_const
p(j) =  log(normpdf(x(j),0.0,0.5)); j=j+1 ;

% chi_p_f
p(j) =  log(betapdf(x(j),3.55,31.75)); j=j+1 ;

% chi_w_f
p(j) =  log(betapdf(x(j),3.55,31.75)); j=j+1 ;

% chi_p
p(j) =  log(betapdf(x(j),3.55,31.75)); j=j+1 ;

% chi_w
p(j) =  log(betapdf(x(j),3.55,31.75)); j=j+1 ;

fprior = sum(p);