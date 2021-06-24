function f = invgam_pdf (x, a, b)
% PURPOSE: returns the pdf at x of the invgamma(a,b) distribution
%---------------------------------------------------
% USAGE: pdf = invgamm_pdf(x,a)
% where: x = vector  
%        a = shape scalar for gamma(a)
%        b = scale 
%---------------------------------------------------
% RETURNS:
%        a vector of pdf at each element of x of the inverse gamma(a,b) distribution      
% --------------------------------------------------
% 
% --------------------------------------------------

%   Code by Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg
%
%   Modified by Mariano Kulish, June-2013

if nargin ~= 3
error('Wrong # of arguments to invgamm_pdf');
end;

if any(any(a<=0))
   error('invgamm_pdf: parameter a is wrong')
end

if any(any(b<=0))
   error('invgamm_pdf: parameter b is wrong')
end

f = (b.^a) .* x .^ (-a-1) .* exp(-b./x) ./ gamma(a);
I0 = find(x<=0);
f(I0) = zeros(size(I0));
