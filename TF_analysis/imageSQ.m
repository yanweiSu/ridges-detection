% imageSQ.m
%
% Display time-frequency result of the Synchrosqueezing transform and others
% function imageSQ(t, ytic, M) ;
%
function imageSQ(t, ytic, M, Qv) ;


S = size(M);
Q = M(:);%Reshape the matrix into a column vector.
Q = Q(find(abs(Q)>0)) ;	

	% truncate the upper bound
q = quantile(Q, Qv);
M(find(M>q)) = q;
%M = M ./ q ;

	% truncate the lower bound
%m = quantile(Q, 0.002) ;
%M(find(M<m)) = m ;


imagesc(t, ytic, M)
axis xy ;
%set(gca, 'fontsize', fz);
