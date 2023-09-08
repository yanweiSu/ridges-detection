function y = ARMA11(N) 
%
% can be changed to ARMA p,q with parameters if you like....
%
dis = random('t', 4, N, 1) ; 
e = armaxfilter_simulate(dis, .0, 1, .5, 1, -.5) ; 
y = e ./ std(e) ;

end