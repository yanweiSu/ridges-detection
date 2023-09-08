function [recon] = Recon_sqSTFT(tfrsq, tfrsqtic, Hz, c, Band)

% Band = absolute frequency band
% tfrsqtic = in the scale of 1Hz

alpha = tfrsqtic(2)-tfrsqtic(1) ;
RR = floor(Band/(Hz*alpha)) ;

recon = [] ;

for kk = 1: length(c)
	idx = max(1,c(kk)-RR): min(length(tfrsqtic),c(kk)+RR) ;
	recon(kk) = 2*sum(tfrsq(idx,kk),1)*(Hz/2/0.5)*(alpha)/Hz ;
end
