// ========================================================================
// Small_economy_dynare.mod
// ========================================================================

// Yield curve definitions
// Define number of yields to include for each country
@# define include_yields = 1
@# define number_of_yields = 8   
@# define observed_yields = [8]

// ======================================================
// MAIN CODE:

// Declare variables

// *** Large economy variables ***

var lambda_f, y_f, infl_f, dy_f, xi_f, z, a_f ;
var n_f, w_f, hp1_f, hp2_f, mc_f, xi_p_f  xi_w_f, hw1_f, hw2_f, winf_f;
var g_f c_f, winf_obs, winf_f_obs ;

varexo eps_xi_f eps_r_f eps_z eps_a_f ;
varexo eps_xi_p_f eps_g_f eps_xi_w_f ;

parameters bet mu rho_r_f phi_pi_f phi_g_f phi_y_f  rho_xi_f rho_z, rho_a_f ;
parameters h_f h psi_f sig_z;
parameters sig_r_f  sig_xi_f, sig_a_f ;
parameters gamm, psi, theta_w_f, eps_w_f, chi_w_f, theta_p_f, eps_p_f, chi_p_f ;
parameters rho_xi_w_f rho_xi_p_f sig_xi_w_f sig_xi_p_f ;
parameters CFoYF rho_g_f sig_g_f ;
// h

// *** Small economy variables ***

var y, dy, de, gamma_F, gamma_X, gamma_H, xi, c g, a;
var infl, pi_F, pi_X, pi_H, lambda, q, rp, r_F ;
var b_F, y_F, c_H, hx1, hx2, x ;
var hf1, hf2, n, mc, w, hp1, hp2,xi_X, xi_H, xi_F  xi_w, hw1, hw2, winf;

varexo  eps_xi eps_r eps_rp eps_g, eps_a ;
varexo eps_xi_X, eps_xi_F, eps_xi_H eps_xi_w ;

parameters alph rho_r phi_pi phi_g phi_y  rho_xi rho_rp sig_a sig_xi sig_rp sig_r, rho_a ;
parameters pit_d pit_f ;
parameters psi_H tau, theta_x, eps_X, chi_X, theta_F, eps_F, chi_F ;
parameters CHoY XoY ;
parameters rho_xi_X, sig_xi_X, rho_xi_F, sig_xi_F, rho_xi_w, sig_xi_w, rho_xi_H, sig_xi_H, sig_a ;
parameters theta_w, eps_w, chi_w, theta_p, eps_p, chi_p ;
parameters dx_const dy_F_const ;
parameters rho_g sig_g ;

// *** Yield curve variables ***

// Declare yields (including policy rates)
@# for j in 1 : number_of_yields
var r_@{j} r_f_@{j};
@# endfor

// Observed yields and shocks 

@# for j in observed_yields
parameters sig_eps_r_@{j} sig_eps_r_f_@{j} r_@{j}_const r_f_@{j}_const;
parameters rho_tp_@{j} rho_tp_f_@{j} ;
var r_@{j}_obs r_f_@{j}_obs tp_@{j} tp_f_@{j};
varexo eps_r_@{j} eps_r_f_@{j};
@# endfor

var r_f_1_obs, dy_f_obs, dc_f_obs, infl_f_obs;
var r_1_obs, dy_obs, dc_obs, infl_obs;//, winf_obs, winf_f_obs ;
var dx_obs dy_F_obs de_obs ;

// Measurement errors
parameters     sig_y_F_shock sig_x_shock;
varexo     y_F_shock x_shock;

// Calibrated parameters - Shared
bet    = 0.9985;
alph   = 0.25;
gamm   = 1;//0.97 ; // Discounting in Euler equation
rho_z  = 0.0 ;
mu     = 1.0043438 ; 

// Calibrated parameters - Large economy 
pit_f		= 1.005 ;
rho_a_f		= 0.0 ;
sig_a_f		= 0.0 ; // Turn off temporary TFP shock
psi_f       = 2 ; 
eps_w_f     = 6 ;
eps_p_f     = 6 ;

// Calibrated parameters - Small economy 
pit_d		= 1.005 ;
rho_a       = 0.0 ;
sig_a		= 0.0 ; // Turn off temporary TFP shock
psi         = 2 ; 
eps_w		= 6 ;
eps_p		= 6 ;
eps_X		= 11 ;
chi_X		= 0 ;
eps_F		= 11 ;
chi_F		= 0 ;
psi_H       = 0.005 ; // Risk premium on NFA position
CHoY        = 0.575 ; 
XoY         = alph ;
CFoYF       = 0.65 ;

// Mode from chain
h_f           =   0.7405 ;
theta_p_f     =   0.7552 ;
theta_w_f     =   0.7819 ;
rho_r_f       =   0.8907 ;
phi_pi_f      =   1.7560 ;
phi_g_f       =   0.0948 ;
phi_y_f       =   0.0817 ;
rho_xi_f      =   0.9512 ;
rho_g_f       =   0.9524 ;
rho_xi_p_f    =   0.9851 ;
rho_xi_w_f    =   0.8033 ;
rho_tp_f_8    =   0.7251 ;
sig_z         =   0.1109 ;
sig_r_f       =   0.1089 ;
sig_xi_f      =   0.2469 ;
sig_g_f       =   0.1263 ;
sig_xi_p_f    =   0.1492 ;
sig_xi_w_f    =   0.1064 ;
sig_eps_r_f_8 =   0.0920 ;
r_f_8_const   =   0.0846 ;
h             =   0.8117 ;
tau           =   3.0631 ;
theta_p       =   0.7948 ;
theta_w       =   0.7792 ;
theta_x       =   0.7774 ;
theta_F       =   0.7663 ;
rho_r         =   0.9112 ;
phi_pi        =   2.2031 ;
phi_g         =   0.1136 ;
phi_y         =   0.1664 ;
rho_rp        =   0.9594 ;
rho_xi        =   0.6377 ;
rho_g         =   0.9294 ;
rho_xi_H      =   0.8488 ;
rho_xi_w      =   0.2944 ;
rho_xi_X      =   0.9091 ;
rho_xi_F      =   0.9620 ;
rho_tp_8      =   0.7388 ;
sig_r         =   0.1664 ;
sig_rp        =   0.3077 ;
sig_g         =   0.1391 ;
sig_xi        =   0.2418 ;
sig_xi_H      =   0.3654 ;
sig_xi_w      =   0.4978 ;
sig_xi_X      =   1.2864 ;
sig_xi_F      =   1.0371 ;
sig_eps_r_8   =   0.1193 ;
r_8_const     =   0.1166 ;
chi_p_f       =   0.0344 ;
chi_w_f       =   0.0808 ;
chi_p		  =   0.0415 ;
chi_w		  =   0.0716 ;

% @# for j in observed_yields
% sig_eps_r_@{j}   = 0.16;
% sig_eps_r_f_@{j} = 0.16;
% 
% r_@{j}_const   = 0.05;
% r_f_@{j}_const = 0.05;
% 
% rho_tp_@{j}      = 0.5 ;
% rho_tp_f_@{j}    = 0.5 ;
% @# endfor

dx_const    = 0 ;
dy_F_const  = 0 ;

sig_y_F_shock     = 0; % 0.88 ;  // Calibrated measurement error for exports
sig_x_shock       = 0; % 0.77 ;  // Calibrated measurement error for imports

//var vd ;

//*******************************//
//*** Model equations ***********//
//*******************************//

model;//(linear) ;

//#h = h_f ; 

// Large economy Taylor rule 9
r_f_1 = rho_r_f*r_f_1(-1) + (1 - rho_r_f)*(phi_pi_f*infl_f +  phi_y_f*y_f) + phi_g_f*dy_f + sig_r_f*eps_r_f/100;

// Small economy, Taylor rule 47
r_1 = rho_r*r_1(-1) + (1 - rho_r)*(phi_pi*infl  + 
    phi_y*y)  + phi_g*dy+ sig_r*eps_r/100;

//*******************************//
//*** Large economy equations ***//
//*******************************//

// Large economy habits 1 - replace y_f with c_f because introduced g
//(mu - bet*h_f)*(mu - h_f)*lambda_f = mu*h_f*y_f(-1) + mu*bet*h_f*y_f(+1)
//    - (mu^2 + bet*h_f^2)*y_f - mu*h_f*(1 - bet*rho_z)*z + (mu - h_f)*(mu - bet*h_f*rho_xi_f)*xi_f;

(mu - bet*h_f)*(mu - h_f)*lambda_f = mu*h_f*c_f(-1) + mu*bet*h_f*c_f(+1)
    - (mu^2 + bet*h_f^2)*c_f - mu*h_f*(1 - bet*rho_z)*z + xi_f;

// Large economy Euler equation with discounting 2
0 = lambda_f - gamm*lambda_f(+1) - ( r_f_1 - infl_f(+1)) + rho_z*z;

// Large economy, hw1_f definition 3
hw1_f = (1-bet*theta_w_f)*(1+psi_f)*(xi_f+n_f) + bet*theta_w_f*(hw1_f(+1) + (1+psi_f)*eps_w_f*(winf_f(+1)-chi_w_f*winf_f)) ;

// Large economy, hw2_f definition 4
hw2_f = (1-bet*theta_w_f)*(lambda_f + w_f + n_f) + bet*theta_w_f*(hw2_f(+1) + (eps_w_f-1)*(winf_f(+1)-chi_w_f*winf_f)) ;

// Large economy, wage index 5
winf_f = (1-theta_w_f)/theta_w_f * (1/(1+eps_w_f*psi_f)) * (hw1_f - hw2_f) + chi_w_f*winf_f(-1)  + xi_w_f;

// Large economy, wage inflation definition 6
winf_f = w_f - w_f(-1) + infl_f + z ;

//(1+psi_f)*(xi_f+n_f) = lambda_f + w_f + n_f ;

// Large economy, aggregate production 7
y_f =  a_f + n_f ; % 

// Large economy, marginal costs 8
mc_f = w_f - a_f ; % 

// Large economy, hp1_f definition 9
hp1_f = (1-bet*theta_p_f)*(lambda_f + mc_f + y_f) + bet*theta_p_f*(hp1_f(+1) + eps_p_f*(infl_f(+1) - chi_p_f*infl_f)) ;

// Large economy, hp2_f definition 10
hp2_f = (1-bet*theta_p_f)*(lambda_f + y_f) + bet*theta_p_f*(hp2_f(+1) + (eps_p_f-1)*(infl_f(+1) - chi_p_f*infl_f)) ;

// Large economy, price index 11
infl_f = ((1-theta_p_f)/theta_p_f)*(hp1_f - hp2_f ) + chi_p_f*infl_f(-1) + xi_p_f ;

// Large economy output growth 10
dy_f = y_f - y_f(-1) + z;

// Large economy demand shock 11
xi_f = rho_xi_f*xi_f(-1) + sig_xi_f*eps_xi_f/100;

// Large economy temp tech shock 12
a_f = rho_a_f*a_f(-1) + sig_a_f*eps_a_f/1000;

// Permanent tech shock 13
z = rho_z*z(-1) + sig_z*eps_z/100;

// Large economy wage shock 17
xi_w_f = rho_xi_w_f*xi_w_f(-1) + sig_xi_w_f*eps_xi_w_f/100;

// Large economy price shock 14
xi_p_f = rho_xi_p_f*xi_p_f(-1) + sig_xi_p_f*eps_xi_p_f/100;

// Large economy government spending shock 15
g_f = rho_g_f * g_f(-1) + sig_g_f * eps_g_f / 10 ;

// Large economy market clearing 16
y_f = CFoYF * c_f + (1 - CFoYF) * g_f ;

//*******************************//
//*** Small economy equations ***//
//*******************************//

// Small economy, habits 21
//(mu - bet*h)*(mu - h)*lambda = mu*h*c(-1) + mu*bet*h*c(+1) - (mu^2 + bet*h^2)*c
//        - mu*h*(1 - bet*rho_z)*z + (mu - h)*(mu - bet*h *rho_xi)*(xi+0*xi_f);

(mu - bet*h)*(mu - h)*lambda = mu*h*c(-1) + mu*bet*h*c(+1) - (mu^2 + bet*h^2)*c
        - mu*h*(1 - bet*rho_z)*z + xi+0*xi_f;

// Small economy, Euler equation with discounting 22
0 = lambda - gamm*lambda(+1) - (r_1 - infl(+1)) + rho_z*z;

// UIP 23
r_1 = r_F + de(+1) ;

// Real exchange rate 24
q = q(-1) + de + infl_f - infl;

// Interest on foreign borrowing 25
r_F = r_f_1 - psi_H * b_F + rp ;

// Imports demand 26
y_F = (c + g) - tau*gamma_F ;

// Demand for home-produced goods 27
c_H = c - tau*gamma_H ;

// Home consumer price index 28
infl = (1 - alph)*(pi_H + gamma_H(-1)) + alph*(pi_F + gamma_F(-1));

// Small economy, definition of hx1 29
hx1 = (1-bet*theta_x)*(lambda + gamma_H + x) + bet*theta_x*(hx1(+1) + eps_X*(pi_X(+1) + chi_X * pi_X)) ;

// Small economy, definition of hx2 30
hx2 = (1-bet*theta_x)*(lambda + gamma_X + x) + bet*theta_x*(hx2(+1) + (eps_X-1)*(pi_X(+1) + chi_X * pi_X)) ;

// Exports inflation 31
pi_X = (1-theta_x)/theta_x * (hx1 - hx2 ) + chi_X * pi_X(-1) + xi_X ;

// Exports demand 32
x = y_f - tau * (gamma_X - q) ;

// Small economy, definition of hf1 33
hf1 = (1-bet*theta_F)*(q + lambda + y_F) + bet*theta_F*(hf1(+1) + eps_F*(pi_F(+1) - chi_F * pi_F)) ;

// Small economy, definition of hf2 34
hf2 = (1-bet*theta_F)*(gamma_F + lambda + y_F) + bet*theta_F*(hf2(+1) + (eps_F-1)*(pi_F(+1) - chi_F * pi_F)) ;

// Imports inflation 35
pi_F = (1-theta_F)/theta_F * (hf1 - hf2 ) + chi_F * pi_F(-1) + xi_F;

// Small economy, production of domestic goods 36
y =  a + n ; 

// Small economy, marginal costs 37
mc = w - gamma_H - a ; % 

// Small economy, definition of hp1 38
hp1 = (1-bet*theta_p)*(lambda + mc + gamma_H + y) + bet*theta_p*(hp1(+1) + eps_p*(pi_H(+1) - chi_p*pi_H)) ;

// Small economy, hp2 definition 39
hp2 = (1-bet*theta_p)*(lambda + gamma_H + y) + bet*theta_p*(hp2(+1) + (eps_p-1)*(pi_H(+1) - chi_p*pi_H)) ;

// Small economy, price index 40
pi_H = ((1-theta_p)/theta_p)*(hp1 - hp2 ) + chi_p*pi_H(-1) + xi_H;

// Small economy, hw1 definition 41
hw1 = (1-bet*theta_w)*(1+psi)*(n+xi+0*xi_f) + bet*theta_w*(hw1(+1) + (1+psi)*eps_w*(winf(+1)-chi_w*winf)) ;

// Small economy, hw2 definition 42
hw2 = (1-bet*theta_w)*(lambda + w + n) + bet*theta_w*(hw2(+1) + (eps_w-1)*(winf(+1)-chi_w*winf)) ;

// Small economy, wage index 43
winf = (1-theta_w)/theta_w * (1/(1+eps_w*psi)) * (hw1 - hw2 ) + chi_w*winf(-1) + xi_w;

// Small economy, wage inflation definition 44
winf = w - w(-1) + infl + z ;

//(1+psi)*(n+xi) = lambda + w + n ;

// Small economy, domestic market clearing 45 - adjusted for g
y = CHoY*c_H + XoY*x + (1 - CHoY - XoY) * (g - tau * gamma_H);

// Small economy, financial market clearing 46
b_F = b_F(-1) / bet + alph * (gamma_X + x - q - y_F) ;

// Exports relative price 48
gamma_X = gamma_X(-1) + pi_X - infl + de;

// Imports relative price 49
gamma_F = gamma_F(-1) + pi_F - infl;

// Home goods relative price 50
gamma_H = gamma_H(-1) + pi_H - infl;

// Home output growth 51
dy = y - y(-1) + z;

// Home temp tech 50
a = rho_a*a(-1) + sig_a*eps_a/1000;

// Home demand 51
xi = rho_xi*xi(-1) + sig_xi*eps_xi/100;

// Risk premium 52
rp = rho_rp * rp(-1) + sig_rp*eps_rp/100;

// Small economy export price shock 53
xi_X = rho_xi_X*xi_X(-1) + sig_xi_X*eps_xi_X/100;

// Small economy import price shock 54
xi_F = rho_xi_F*xi_F(-1) + sig_xi_F*eps_xi_F/100;

// Small economy domestic price shock 55
xi_H = rho_xi_H*xi_H(-1) + sig_xi_H*eps_xi_H/100;

// Small economy wage shock 56
xi_w = rho_xi_w*xi_w(-1) + sig_xi_w*eps_xi_w/100;

// Small economy government spending shock
g = rho_g * g(-1) + sig_g * eps_g / 10 ; 

//****************************//
//*** Bond yield equations ***//
//****************************//

// Additional bond yields
@# if include_yields == 1
@# for j in 2 : number_of_yields
r_f_@{j} = (1/@{j})*(r_f_1 + (@{j}-1)*r_f_@{j-1}(+1));
r_@{j}   = (1/@{j})*(r_1 + (@{j}-1)*r_@{j-1}(+1));
@# endfor
@# endif

//****************************//
//*** Observable equations ***//
//****************************//

// Observation equations for yields
@# for j in observed_yields
r_@{j}_obs   = 400 * (log(pit_d) + log(mu) - log(bet) + r_@{j} + r_@{j}_const/100 + tp_@{j});
r_f_@{j}_obs = 400 * (log(pit_f) + log(mu) - log(bet) + r_f_@{j} + r_f_@{j}_const/100 + tp_f_@{j});

// Term premia
tp_@{j}     = rho_tp_@{j} * tp_@{j}(-1) + sig_eps_r_@{j} * eps_r_@{j}/100;
tp_f_@{j}   = rho_tp_f_@{j} * tp_f_@{j}(-1) + sig_eps_r_f_@{j} * eps_r_f_@{j}/100;


@# endfor

r_1_obs   = 400 * (r_1 + (log(pit_d) + log(mu) - log(bet))) ;
r_f_1_obs = 400 * (r_f_1 + (log(pit_f) + log(mu) - log(bet))) ;

dy_obs = 100 * (dy + log(mu)) ;
dy_f_obs = 100 * (dy_f + log(mu)) ;

infl_obs = 100 * (infl + log(pit_d))  ;
infl_f_obs = 100 * (infl_f + log(pit_f))  ;

winf_obs = 100 * (w - w(-1) + log(pit_d) + log(mu)) ;
winf_f_obs = 100 * (w_f - w_f(-1) + log(pit_f) + log(mu));

dx_obs = 100 * (x - x(-1) + z + log(mu) + dx_const) + sig_x_shock * x_shock  ;
dy_F_obs = 100 * (y_F - y_F(-1) + z + log(mu) + dy_F_const) + sig_y_F_shock * y_F_shock  ;

de_obs = 100 * (de + log(pit_d) - log(pit_f)) ; 

dc_obs   = 100 * (c - c(-1) + z + log(mu)) ; 
dc_f_obs = 100 * (c_f - c_f(-1) + z + log(mu)) ; 

//vd = exp(xi)*(log(exp(c)-h*exp(c(-1))) - 1/(1+psi) * exp(n)^(1+psi)) + bet*vd(+1) ;  

end;

% // *** Initialise variables

initval;

lambda_f    = 0;
y_f         = 0;
infl_f      = 0;
z           = 0;
xi_f        = 0;
lambda      = 0;
c           = 0;
xi          = 0;
q           = 0;

end;

// Initialise shocks
shocks;

var eps_xi_f;   stderr 1;
var eps_a_f;    stderr 1;
var eps_xi_w_f; stderr 1;
var eps_xi_p_f; stderr 1;
var eps_r_f;    stderr 1;
var eps_g_f;    stderr 1;

var eps_xi;     stderr 1;
var eps_a;      stderr 1;
var eps_rp;     stderr 1;
var eps_xi_X;   stderr 1;
var eps_xi_F;   stderr 1;
var eps_xi_w;   stderr 1;
var eps_xi_H;   stderr 1;
var eps_r;      stderr 1;
var eps_g;      stderr 1;

@# for j in observed_yields
var eps_r_@{j}; stderr 1;
var eps_r_f_@{j}; stderr 1;
@# endfor

var eps_z ; stderr 1 ;

//var n_shock ; stderr 1 ;
//var n_f_shock ; stderr 1 ;

var y_F_shock ;     stderr 1 ;
var x_shock ;       stderr 1 ;

end;

steady;

stoch_simul(order=1,nograph) r_f_1_obs r_f_8_obs dy_f_obs dc_f_obs infl_f_obs winf_f_obs r_1_obs r_8_obs dy_obs dc_obs infl_obs winf_obs dy_F_obs dx_obs de_obs ;
//stoch_simul(order=1,irf=20) dy_obs dy_f_obs r_1_obs r_f_1_obs r_8_obs r_f_8_obs infl_obs infl_f_obs winf_obs winf_f_obs dx_obs dy_F_obs dc_obs dc_f_obs de_obs;
// stoch_simul(order=1,periods=120,nograph);

