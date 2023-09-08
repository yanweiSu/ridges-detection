clear;
close all;
addpath('TF_analysis/');
addpath(genpath('./MFEToolbox/'));
addpath('./RRP-RD/toolbox/TF_Toolbox/');
addpath('./RRP-RD/RRP_alg/');

%%
AFUND2 = 1.0;
SNR = [Inf 5 0];
AFUND1 = [.1 0.2 0.5];

load ./newARMAnoise/curve_MultiCurveExt.mat
load ./newARMAnoise/RMSE_MultiCurveExt.mat
load ./newARMAnoise/time_MultiCurveExt.mat

% curve_MultiCurveExt = {};
% RMSE_MultiCurveExt = {};
% time_MultiCurveExt = {};
% 
% curve_RRPRD = {};
% RMSE_RRPRD = {};
% time_RRPRD = {};

for am = 2:-1:1

startSNR = 1;
if am == 3
    startSNR = 3;
end

for snrdb = startSNR:3
% am = 1; snrdb = 3; r = 76;

% Generate signal sample
start = 1;
curveMH = zeros(2048,100);
rmseMH = zeros(100,1);
elapsedMH = zeros(100,1);

if snrdb == 3 && am == 3
    start = 66;
    curveMH = curve_MultiCurveExt{am,snrdb};
    rmseMH = RMSE_MultiCurveExt{am,snrdb};
    elapsedMH = time_MultiCurveExt{am,snrdb};
end

D = AFUND1(am);
snr = SNR(snrdb);

for r = start:100
fprintf("D = %.1f, SNR = %d, r = %d; ", D, snr, r);

% change this random seed if you want to see other results
rand('state', r);
randn('state', r);
%% generate simulated signals 

% start with a high sampling rate so that the phase can be more
% accurately generated (will downsample it later)
Hz = 2000 ;
t = (1:1e5)' / Hz ;

% generate AM for the 1st component
AM1 = smoothdata(abs(cumsum(randn(size(t)))), "gaussian", .4*length(t));
AM1 = 3*AM1 ./ max(AM1) + 2.5;
AM1 = AM1 .* exp(-((t-10)/30).^2);

% generate FM for the 1st component
IF1r = smoothdata(abs(cumsum(randn(size(t)))), "gaussian", .4*length(t));
% IF1 is between pi/2+0.4+/-1; that is, around [1,3] Hz
IF1r = IF1r./max(abs(IF1r)) + .97;
% IF1 = smooth(IF1r, 8000);

% this is the smoothing step. Can use the shorter window to get more complicated IF
tmp = cumsum(exp(-((t-25)/2.5).^2));
tmp = 1.5*tmp ./ max(tmp);
IF1 = IF1r + tmp;

phi1 = [];
for l = 1:3
    u = randn(size(t))*sqrt(0.1);
    u = smoothdata(u, "gaussian", 0.1*length(t));
    phi1 = [phi1 l*(cumsum(IF1)./Hz + t.^1.9/(20*1.9)) + u];
end

% generate the 1st oscillatory component (like PPG)
aa = rand(1,2); %a4 = aa(1);
a3 = aa(1);
a2 = sum(aa);
M = max([D a2 a3]);
a1 = D./M; a2 = a2./M; a3 = a3./M; %a4 = a4./M;
s3 = a1 * cos(2*pi*phi1(:,1)) ...
    + a2 * cos(2*pi*phi1(:,2)) ...
    + a3 * cos(2*pi*phi1(:,3));
%     + a4 * cos(2*pi*4*phi1+ph(4));
x1 = AM1 .* s3;

% generate AM & FM for the 2nd component
AM2 = smooth(abs(cumsum(randn(size(t)))./Hz), 6000) ;
AM2 = 3*AM2 ./ max(AM2) + 2.3 ;
AM2 = AM2 .* exp(-((-8e4+1:2e4)'/5e4).^1.8) ;
ff = abs(cumsum(randn(1e5, 1))) ;   % standard brownian motion
% IF2 is shifted up from IF1 by exp(1)/2 +/- 1 ~ [0.35, 2.35]
IF2r = IF1r + ff./max(abs(ff)) + 1.36 ;
% this is the smoothing step. Can use the shorter window to get more
% complicated IF
IF2 = smooth(IF2r, 5000) + t/20 ;   % smooth(standard brownian) + e/2 + t/20
phi2 = cumsum(IF2) ./ Hz + t.^2/(20*2) ;


%% generate the 2nd oscillatory component (with the very general shape) you can change the shape to contain only 5 or 8 harmonics
% gg = mod(phi2,1);
% [a,b] = findpeaks(gg);
% b = [1; b; 2*b(end)-b(end-1)] ;
% s2 = zeros(size(phi2)) ;
% for ii = 1: length(b)-1
%     idx = b(ii):b(ii+1) ;
%     s2(idx) = (idx-b(ii)) ./ (b(ii+1)-b(ii)+1) ;
% end
% x2 = AM2 .* s2(1:length(AM2)) ;
% x2 = x2 - mean(x2)/2 ;
% 
% % generate the 2nd oscillatory component (like PPG)
% bb = rand(1,4); b4 = bb(4); b3 = sum(bb(3:4)); b2 = sum(bb(2:4)); b1 = sum(bb);
% M = max([b1 b2 b3 b4]); b1 = AFUND2*b1./M ; b2 = AFUND2*b2./M ; b3 = AFUND2*b3./M;
% ph = rand(1,4)*2*pi-pi;
% s3 = b1.*cos(2*pi*phi2+ph(1)) + b2.*cos(2*pi*2*phi2+ph(2)) ...
%     + b3.*cos(2*pi*3*phi2+ph(3)) + 0.6*b3.*cos(2*pi*4*phi2+ph(4)); 
% x2 = AM2 .* s3;

%%
x = x1; % + x2;

% % Old noise
% noise = random('T',5,length(x),1) ;
% sigma = 10^(log10(std(x)./std(noise)) - snr/20);
% noise = sigma * noise ; 
% % var(noise)
% % snrdb = 20 * log10(std(x)./std(noise)) ;
% fprintf(['SNR = ',num2str(20*log10(std(x)./std(noise))),'\n']);

% generate 600 points (depends on your situation) for the first part.
noise = zeros(length(x),1);
% y = armaxfilter_simulate(length(x)/2, .0, 2, [0 -.95], 1, -.5);
y = ARMA11(length(x)/2);
% sigma = 10^(log10(std(x(1:length(x)/2))./std(y)) - (snr)/20);
noise(1:length(x)/2) = y;   % * sigma;

y = random('T', 5, length(x)/2, 1);
% sigma = 10^(log10(std(x(length(x)/2+1:end))./std(y)) - (snr)/20);
noise(length(x)/2+1:end) = y;   % * sigma;

sigma = 10^(log10(std(x)./std(noise)) - (snr)/20);
noise = noise * sigma;

fprintf("SNR = %.1f db\n", 20*log10(std(x)./std(noise)));

Y = x + noise;

% downsample to 200Hz.
x = x(1:10:end) ;
x1 = x1(1:10:end);
% x2 = x2(1:10:end);
Y = Y(1:10:end) ;
AM1 = AM1(1:10:end) ;
IF1 = IF1(1:10:end) ;
AM2 = AM2(1:10:end) ;
IF2 = IF2(1:10:end) ;
t = t(1:10:end) ;
noise = noise(1:10:end) ;
orgHz = Hz;
Hz = Hz / 10;

fs = 50;
x = resample(x, fs, 200);  % resample again to process faster
Y = resample(Y, fs, 200);
x1 = resample(x1, fs, 200);
% x2 = resample(x2, 50, 200);

% IF = 2000*diff(phi2);
% IF = resample(IF, fs, 2000);
IF1 = 2000*diff(phi1(:,1));
IF1 = resample(IF1, fs, 2000);
for k = 2:3
    tmp = 2000*diff(phi1(:,k));
    tmp = resample(tmp, fs, 2000);
    IF1 = [IF1 tmp];
end

%%
len = 2^(floor(log2(length(Y))));
xn = Y(1:len);

IF1_fund = IF1(1:len,1);
% c_est = curve_RRPRD{am,snrdb}(:,r);
% rmse(r) = norm((IF1_fund-c_est),2)/norm(IF1_fund,2);
% 
% end
% 
% RMSE_RRPRD{am,snrdb} = rmse;
% save('./RMSE_RRPRD.mat', 'RMSE_RRPRD');
% end
% end

% figure; plot(IF1_fund); hold on; plot(curve_RRPRD{am,snrdb}(:,r)); hold off;

%% signal samples
Lx = length(xn);
Tx = (0:length(xn)-1)./fs;

%% time frequency parameters
sigma_opt = 1/sqrt(600); % minimizing renyi entropy
[g, ~] = gauss_win(Lx, sigma_opt);

% second order synchrosqueezing & ridge detection
Nfft = length(xn)/2;
[~,TFR] = sst2(xn, sigma_opt, Nfft);
fr = fs/Nfft;

% tfrsq for MHRD
tfrsq = TFR.SST2(1:Nfft/2,:);
Fx = (0:Nfft/2-1)*fr;

% % resample for MHRD
% tfr_res = tfrsq(1:2:end,:);
% fr_res = 2*fr;
% Fx_res = (0:Nfft/4) * fr_res;

% figure;
% imageSQ(Tx, Fx, abs(tfrsq), 0.99); ylim([0 15])
% colormap(1-gray); colorbar;
% set(gca,'YDir','normal');
% xlabel("time");
% ylabel("frequency");
% title("short time Fourier transform");
% 
% figure; plot(real(Y)); hold on; plot(real(x)); hold off;

% figure;
% imageSQ(Tx, Fx_res, abs(tfr_res), 0.99); ylim([0 15])
% colormap(1-gray); colorbar;
% set(gca,'YDir','normal');
% xlabel("time");
% ylabel("frequency");
% title("resample");


%% SST by Marcelo
% [STFT, SST1, SST2, ~, ~, ~, ~, ~, tau2, tau3, phi22p, ~, ~] ...
%     = sstn_test_new(xn, 0, sigma_opt, Nfft, 1:length(xn));
% fr = fs/Nfft;
% Fx = (0:Nfft/2-1)*fr;
% figure;
% imageSQ(Tx, Fx, abs(SST2), 0.99); ylim([0 15])
% colormap(1-gray); colorbar;
% set(gca,'YDir','normal');
% xlabel("time");
% ylabel("frequency");
% title("short time Fourier transform");

%% RD
tfr = abs(tfrsq);
tfrtic = Fx;
tfrtic = reshape(tfrtic, [], 1);

lambda = [2.0 1.8 1.6];
multi = [2 3];
mu = zeros(size(multi));
band_fund = [.3 4.5];
band_multi = [.2 .3];
% Inital: larger searching band and longer time
tt = [5*fs 10*fs:fs:35*fs];
% tt = 0:fs:size(tfr,2); tt = [tt size(tfr,2)];

cALL = [];
flag = 0;
fprintf("Run 3-curves simultaneous extraction\n")
tic
c_init = [1; 1; 1];
for l = 1:length(tt)-1
    TFRtmp = abs(tfr(:, tt(l)+1:tt(l+1)));
    c0 = ones(size(TFRtmp,2),1);
    c1 = ones(size(c0));
    c2 = ones(size(c0));
    fprintf("[%.1fsec, %.1fsec]x[%.2fHz, %.2fHz]:\n", tt(l)/fs, tt(l+1)/fs, band_fund(1), band_fund(2));
    [c0, c1, c2] = CurveMultiExt_init_3curves(TFRtmp.', tfrtic, lambda, mu, multi, band_fund, band_multi, flag, c_init);
    c_init = [c0(end); c1(end); c2(end)];
    band_fund = [max(tfrtic(c0(end))-.4, 0), min(tfrtic(c0(end))+1.0, tfrtic(end))];
    cALL = [cALL; [c0 c1 c2]];
    flag = 1;
end
elapsedMH(r) = toc;

% cALL = MultiCurveExt(tfrsq, tfrtic, lambda, mu, multi, band_fund, band_multi, 2, tt, [-.5, 1.0]);

% padding ones
cALL_pad = ones(length(xn), size(cALL,2));
cALL_pad(tt(1)+1:tt(end), :) = cALL;
% figure; plot(IF1_fund); hold on; plot(Fx_res(cALL_pad(:,1))); hold off;

curveMH(:,r) = tfrtic(cALL_pad(:,1));
seg = tt(1)+1 : tt(end);
rmseMH(r) = norm((tfrtic(cALL_pad(seg,1))-IF1_fund(seg)),2)/norm(IF1_fund(seg),2);
fprintf("%d, RMSE = %.2f, time = %.2f\n", r, rmseMH(r), elapsedMH(r));

curve_MultiCurveExt{am,snrdb} = curveMH;
RMSE_MultiCurveExt{am,snrdb} = rmseMH;
time_MultiCurveExt{am,snrdb} = elapsedMH;

save('./newARMAnoise/curve_MultiCurveExt.mat', 'curve_MultiCurveExt');
save('./newARMAnoise/RMSE_MultiCurveExt.mat', 'RMSE_MultiCurveExt');
save('./newARMAnoise/time_MultiCurveExt.mat', 'time_MultiCurveExt');

end
end
end

%%
% for t = 1:length(tt)-1
%     TFR = tfr(:, tt(t)+1:tt(t+1)); %CurveMultiExt_3curves
%     tic
%     fprintf("%.0f to %.0f sec, %.2f Hz ~ %.2f Hz:\n", tt(t)/fs, tt(t+1)/fs, band_fund(1), band_fund(2));
%     [c0, c1] = CurveMultiExt_init_2curves_new(TFR.', tfrtic.', ...
%         lambda, mu, multi, band_fund, band_multi, flag, c_init);
%     toc
%     c_init = [c0(end); c1(end)];
%     band_fund = [max(tfrtic(c0(end))-1.0, 0), min(tfrtic(c0(end))+1.0, fs/2)];
% %     band_fund = [tfrtic(c0(end))*fs - 0.3, tfrtic(c0(end))*fs + 0.5];
%     cALL = [cALL; [c0 c1]];
%     flag = 1;
% end

%%
% load curve_RRPRD.mat
% figure; plot(IF1_fund); hold on;
% plot(curve_RRPRD{am,snrdb}(:,r)); hold on;
% plot(Fx(cALL(:,1))); hold off;
% 
% seg = 10*fs+1:40*fs;
% norm(curve_RRPRD{am,snrdb}(seg,r)-IF1_fund(seg))/norm(IF1_fund(seg))
% norm(Fx(cALL(seg,1)).'-IF1_fund(seg))/norm(IF1_fund(seg))
% 
% 
% toc
% end

% % save('./curve_RRPRD.mat', 'curve_RRPRD');
% end
% end

%% 
% figure;
% imageSQ(Tx, Fx_res, abs(tfr_res), 0.99); ylim([0 15]);
% colormap(1-gray); colorbar;
% set(gca,'YDir','normal');
% xlabel("time"); ylabel("frequency");
% for r = 1:3
%     hold on;
%     plot(Tx, Fx_res(cALL_pad(:,r)));
% end
% hold off;
% title("result of the ridge detection");

% comb = abs(tfrsq);
% for t = 1:size(comb,2)
% for f = 1:size(comb,1)
% for k = 2:3
%     RR = max(1, k*f-10) : min(size(comb,1), k*f+10);
%     add = max(abs(tfrsq(RR,t)));
%     if isempty(RR)
%         add = 0;
%     end
%     comb(f,t) = comb(f,t) + add;
% end
% end
% end
% 
% figure;
% imageSQ(Tx, Fx, comb, 0.99); ylim([0 15]);
% colormap(1-gray); colorbar;
% set(gca,'YDir','normal');
% xlabel("time"); ylabel("frequency");
% title("result of the ridge detection");
% 
% idx = find(tfrtic > 0.5 & tfrtic < 5);
% [c] = CurveExt(comb(idx,:).', 1.0);
% c = c + idx(1) - 1;
% hold on;
% plot((0:length(c)-1)./fs, tfrtic(c), 'r');
% hold off;
