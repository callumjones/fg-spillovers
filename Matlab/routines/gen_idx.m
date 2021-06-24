function out = gen_idx(in,t_in)

if nargin<2 ; t_in = 1:length(in) ; end

out(1) = 1 ;

for tt=2:length(in)
    out(tt) = out(tt-1)*(1+in(tt)) ;
end

if length(t_in)>length(in) 
    for tt=1:(find(t_in==1)-1)
        out = [out(1)/(1+in(1)), out] ;
    end
    out = out / out(find(t_in==0)) ;
end