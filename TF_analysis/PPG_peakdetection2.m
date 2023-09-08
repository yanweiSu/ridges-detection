function [R, Rs, Q, M] = PPG_peakdetection2(ppg, Fs)
%PPG Find the R peaks from the signal ppg sampled at Fs Hz.

%%
tic
ppg = ppg(:);
ppg(isnan(ppg)) = 0;
switch nargin
    case 1
        error('Specify the sampling rate Fs as the second argument.')
end

%%
% parameters: 設定thershold; 最小寬度
threshold = 0.75; % 0.02
QRS_length = fix(0.11 * Fs);
beat_length = fix(0.66 * Fs);

% QRS isolation, 用bandpass fillter , v1 = indicatot, v2 = V, mu, z,  V1>v2
% for at least QRS length, 
[b, a] = butter(2, [0.5, 8] / (Fs / 2));
filtered = filtfilt(b, a, ppg);
filtered(filtered < 0) = 0;
u = filtered.^2;
V = filtfilt(ones(1, beat_length) ./ beat_length, 1, u);
indicator = filtfilt(ones(1, QRS_length) ./ QRS_length, 1, u);

mu = smooth(u, Fs*30, 'loess') ;

try threshold = filtfilt(ones(1, 5 * Fs) ./ (5 * Fs), 1, u) * threshold;
catch
    threshold = threshold .* mu ; %* mean(u);
end

t = indicator > V + threshold;

% QRS detection 在M那段裡面，找到local maximum
[M, start] = regexp(sprintf('%i', [0 t']), '1+', 'match'); 
M = cellfun(@length, M);
R = nan(length(M), 1);
for i = 1:length(M)
    if M(i) > QRS_length
        [~, R(i)] = max(ppg(start(i) + 1:min(start(i) + M(i) - 1, length(ppg))));
        R(i) = R(i) + start(i);
    end
end
R0 = R(~isnan(R));

R = R0 ; 

idx = [];
for jj = length(R):-1:2
   if R(jj)-R(jj-1)< ceil(0.4*Fs)
       idx = [idx jj] ;
       
       if jj>2
           if R(jj) - R(jj-2)< ceil(0.4*Fs) 
                idx = [idx jj-1] ;
           end
       end
   end
end

R(unique(idx)) = [] ;

%% get Q points
Q = zeros(size(R)) ;
for jj = 1: length(R)
    idx = (max(1, R(jj)-0.5*Fs):R(jj));
    [~, tmp] = min(ppg(idx)) ;
    %Q(jj) = max(1, R(jj) - 0.2*Fs) + 1 + tmp(1) ;
    Q(jj) = max(1, R(jj)-0.5*Fs) + tmp(1) - 1 ; 
    % Should be -1 so that if tmp = 1 we can still take the correct one.
end

%% get the final peaks (largest slope over the upstroke/systolic phase)
[b, a] = butter(2, [0.5, 12] / (Fs / 2));
filtered = filtfilt(b, a, ppg);
dPPG = derivative([1:length(filtered)]'/Fs, filtered) ;
Rs = zeros(size(R)) ;
for jj = 1: length(R)
    idx = [Q(jj): R(jj)] ;
    [~, tmp] = max(dPPG(idx)) ;
    Rs(jj) = Q(jj) - 1 + tmp(1) ;
end


disp(['Peak detection completed in ' num2str(toc) ' seconds.'])


end


