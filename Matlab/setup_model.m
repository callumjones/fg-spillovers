
model_name = 'Small_economy_dynare' ;

%% Structure 1: high inflation - high debt

eval(['dynare ', model_name, '.mod noclearall nostrict']) ;

dyn_in.M_  = M_ ;
dyn_in.oo_ = oo_ ;
dyn_in.solve = 1 ;
dyn_in.linearize_around_diff_y = 0 ;
dyn_in.options_ = options_ ;

SET.M_ = dyn_in.M_ ;
SET.oo_ = dyn_in.oo_ ;
SET.options_ = options_ ;

out = dyn_to_str(dyn_in) ;

SET.mats.Q = out.mats.Q ; Q = SET.mats.Q ;
SET.mats.G = out.mats.G ; G = SET.mats.G ;

SET.str_mats.A0 = out.mats.A ;
SET.str_mats.A1 = out.mats.B ;
SET.str_mats.B0 = out.mats.D ;
SET.str_mats.D0 = out.mats.E ;
SET.str_mats.D2 = out.mats.D2 ;
SET.str_mats.Gamma = out.mats.Gamma ;

gen_set_params(SET,dyn_in)

%% Extract variable details from Dynare output

n_ = M_.endo_nbr ; % Number of endogenous variables
l_ = M_.exo_nbr ;  % Number of exogenous variables
n_ = n_+1 ;        % Constant

SET.variable.n_ = n_ ; 
SET.variable.l_ = l_ ; 
SET.nparam      = M_.param_nbr ;

var_names   = cellstr(M_.endo_names) ;
param_names = cellstr(M_.param_names) ;
exo_var_names = cellstr(M_.exo_names) ;

for ii=1:n_-1
    eval(['SET.variable.',var_names{ii},'=',int2str(ii),';']) ;
end

for ii=1:l_
    eval(['SET.variable.shock.',exo_var_names{ii},'=',int2str(ii),';']) ;
end

for ii=1:M_.param_nbr
    eval(['SET.param_names.',param_names{ii},'=',int2str(ii),';']) ;
end

SET.param_names_num = SET.param_names ; % problem with cell and strings

%% Compute steady-state

SET.y_ss = (eye(n_-1) - Q(1:n_-1,1:n_-1)) \ Q(1:n_-1,end)  ;
SET.y_ss = [SET.y_ss ; 1] ; % Constant

xsteady_model = SET.y_ss(1:end-1) ;
xsteady=oo_.steady_state ;

% Check steady state from equations and Dynare

fprintf('\n')
fprintf('\n')
fprintf('Check steady-state... sum of SS and model SS is:   %5.4f', sum(xsteady-xsteady_model)) ;
fprintf('\n')
fprintf('\n')
fprintf('This number should be zero')
fprintf('\n')
fprintf('\n')

%% Tidy Up Dynare Files
% 
% %eval(['!mv ' model_name '*.m ./dynareoutput/']) ;
% movefile([model_name, '*.m'],'./dynareoutput/')

eval(['delete ' model_name '*.m']);
eval(['delete ' model_name '*.mat']);
eval(['delete ' model_name '*.swp']);
eval(['delete ' model_name '*.log']);
eval(['delete ' model_name '*.eps']);
eval(['delete ' model_name '*.pdf']);
eval(['delete ' model_name '*.fig']);
eval(['delete ' model_name '*.jnl']);
eval(['delete ' model_name '*.log']);
eval(['rmdir(''./' model_name ''',''s'')'])

delete g*.mat H.dat

%% ZLB algorithm contents

SET.maxchk   = 200 ; % Check ZLB binding out 30 periods
SET.stoch    = 1 ;   % Shocks each period? For IRF, identical to 0

% Rows of foreign and domestic Taylor rules
SET.f_tr = 1 ;
SET.d_tr = 2 ;

SET.i_d = SET.variable.r_1 ;
SET.i_f = SET.variable.r_f_1 ;

SET.param_zlb_f     = log(pit_f) + log(mu) - log(bet) - 0.00125/4 ;
SET.param_zlb_d     = log(pit_d) + log(mu) - log(bet) - 0.0025/4;
SET.param_zlb_d_1pc = log(pit_d) + log(mu) - log(bet) - 0.01/4 ;

mats = resolve(SET,dyn_in) ; % Construct structural matrices: A0, A1, etc

SET.mat_init     = mats.mat_init ;
SET.mat_i_f_zlb  = mats.mat_i_f_zlb ;
SET.mat_i_d_zlb  = mats.mat_i_d_zlb ;
SET.mat_i_df_zlb = mats.mat_i_df_zlb ;
SET.mat_fin      = mats.mat_fin ;
