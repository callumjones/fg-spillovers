function [mean, std] = prior_m2d(a, b, dist)

x0 = [0.05, 0.2];

options=optimset('Display','iter'); 

switch dist
    
    case 'unif'

        mean = 0.5 * (a+b) ;
        std  = sqrt(1/12 * (b-a)^2) ;

    case 'beta'
        
        mean = a / (a+b) ;
        std  = sqrt((a*b) / ((a+b)^2 * (a+b+1))) ;

    case 'gamma'
        
        mean = a*b ;
        std  = sqrt(a*b^2) ;
        
        %g = @(x) [x(1)*x(2) - mean ; ...
        %    x(1) * x(2)^2 - std^2 ] ;
        %[x,~] = fsolve(g,x0,options) ;
        
        %mean = x(1) ;
        %std  = x(2) ;
        
    case 'norm'
        
        mean = a ;
        std  = b ;
        
    case 'invgamma'
        
        mean = b / (a-1) ;
        std  = sqrt(b^2 / ((a-1)^2*(a-2))) ;

end



