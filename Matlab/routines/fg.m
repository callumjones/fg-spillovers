function [out, zlb_f, z_f_T] = fg(SET, in, e_in, FG)
%function out = fg(SET, in, e_in, FG)
%
% callum.jones@nyu.edu

z_f_i = FG.z_f_i ;
z_f_f = FG.z_f_f ;

z_f_T = z_f_f - z_f_i + 1 ;

maxchk = SET.maxchk ;

%% FG

max_t = z_f_f ;

mats = rfmats(SET, z_f_i, z_f_f) ;

Qt = mats.Qhat ;
Gt = mats.Ghat ;

ychk(:,1) = Qt(:,:,1) * in + Gt(:,:,1) * e_in ;

for tt = 2:max_t
    ychk(:,tt) = Qt(:,:,tt) * ychk(:,tt-1) ;
end

FG.ychk = ychk ; 

[tmp, zlb_f, z_f_T_] = zlb(SET, in, e_in, FG) ; 

if z_f_T_>z_f_T ; z_f_T = z_f_T_ ; end

%%

switch SET.stoch
    case 1 ; out = tmp(:,1) ;
    case 0 ; out = tmp(:,1:maxchk) ; 
end