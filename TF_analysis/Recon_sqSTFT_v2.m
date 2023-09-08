function [recon] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, c, Band, coeff) 

% Band = absolute frequency band usually is chosen to be 0.1Hz
% tfrsqtic = in the scale of 1Hz
% c is the extracted curve.

%coeff=0.1152; % computed in the numerical example 'genFig1Fig2_test_coeff_SST.m'

alpha = tfrsqtic(2)-tfrsqtic(1) ;
RR = round(Band/(Hz*alpha)); 


recon = [] ;

C = 2 * alpha / coeff;
% Is coeff h(0)? (h is window func of STFT)

for kk = 1: length(c)
	idx = max(1,c(kk)-RR) : min(length(tfrsqtic),c(kk)+RR) ;
	recon(kk) = C * sum(tfrsq(idx,kk),1) ;
end

