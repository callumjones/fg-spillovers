% FG_Data.m
%
% This file loads and transforms the data for the Forward Guidance project
% without bond yields.
% =========================================================================

% Load the raw data
data = readtable('FG_Data_Dynare.xls','Sheet','Matlab_Yields', ...
    'Range','A1:AA147');

data_dates = data.date ;

Can_GDP_growth = 100.*(log(data.Can_GDP(2:end)) - log(data.Can_GDP(1:end-1)));
Can_C          = 100.*(log(data.Can_C(2:end)) - log(data.Can_C(1:end-1)));
Can_Infl       = 100.*(log(data.Can_Infl(2:end)) - log(data.Can_Infl(1:end-1)));
Can_R1         = data.Can_R1(2:end);
Can_R8         = data.Can_R8(2:end);
Can_dER        = 100.*(log(data.Can_ER(2:end)) - log(data.Can_ER(1:end-1)));
%Can_N          = data.Can_N(2:end) - data.Can_N(1:end-1);
Can_Wage       = 100.*(log(data.Can_Wage(2:end)) - log(data.Can_Wage(1:end-1)));
Can_M          = 100.*(log(data.Can_M(2:end)) - log(data.Can_M(1:end-1)));
Can_X          = 100.*(log(data.Can_X(2:end)) - log(data.Can_X(1:end-1)));

US_GDP_growth = 100.*(log(data.US_GDP(2:end)) - log(data.US_GDP(1:end-1)));
US_C          = 100.*(log(data.US_C(2:end)) - log(data.US_C(1:end-1)));
US_Infl       = 100.*(log(data.US_Infl(2:end)) - log(data.US_Infl(1:end-1)));
US_R1         = data.US_R1(2:end);
US_R8         = data.US_R8(2:end);
%US_N          = data.US_N(2:end) - data.US_N(1:end-1);
US_Wage       = 100.*(log(data.US_Wage(2:end)) - log(data.US_Wage(1:end-1)));

US_R1(133) = 0.125 ; % Ensure ZLB at 0.125% in US

Zobs = [ ...
    US_R1' ;        % 1
    US_R8' ;        % 2
    US_GDP_growth'; % 3
    US_C';          % 4
    US_Infl';       % 5
    US_Wage';       % 6
    Can_R1' ;       % 7
    Can_R8';        % 8
    Can_GDP_growth';% 9
    Can_C';         % 10
    Can_Infl';      % 11
    Can_M' ;        % 12
    Can_X';         % 13
    Can_dER';       % 14
    Can_Wage'];     % 15
 
Zobs = Zobs(:,35:end) ;    % Select post 1991 data

data_dates = data_dates(36:end); 

meangrowth = mean([Zobs(3,1:68),Zobs(9,1:68)]) ;

% Adjust mean of Canadian wages
Zobs(15,:) = Zobs(15,:) - mean(Zobs(15,1:68)) + meangrowth + mean(Zobs(11,1:68)) ;

% Adjust means of Canadian imports and exports
Zobs(12,:) = Zobs(12,:) - mean(Zobs(12,1:68)) + meangrowth ;
Zobs(13,:) = Zobs(13,:) - mean(Zobs(13,1:68)) + meangrowth ;

% Adjust mean of US wages
Zobs(6,:)  = Zobs(6,:) - mean(Zobs(6,1:68)) + meangrowth + mean(Zobs(5,1:68)) ;

% Adjust mean of Canadian and US GDP 
Zobs(9,:)  = Zobs(9,:) - mean(Zobs(9,1:68)) + meangrowth ;
Zobs(3,:)  = Zobs(3,:) - mean(Zobs(3,1:68)) + meangrowth ;

% Adjust mean of Canadian and US Consumptions
Zobs(10,:) = Zobs(10,:) - mean(Zobs(10,1:68)) + meangrowth ;
Zobs(4,:)  = Zobs(4,:) - mean(Zobs(4,1:68)) + meangrowth ;

% Demean nominal exchange rate
Zobs(14,:) = Zobs(14,:) - mean(Zobs(14,1:68)); 
