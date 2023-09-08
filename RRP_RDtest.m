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

% load ./ARMAnoise/curve_RRPRD.mat
% load ./ARMAnoise/RMSE_RRPRD.mat
% load ./ARMAnoise/time_RRPRD.mat
curve_RRPRD = {};
RMSE_RRPRD = {};
elapsed_time = {};

for am = 1:3
for snrdb = 3:-1:1
% am = 1; snrdb = 2;

D = AFUND1(am);
snr = SNR(snrdb);

% Generate signal sample
curves = zeros(2048,100);
rmse = zeros(100,1);
elapsed = zeros(100,1);
% curves = curve_RRPRD{am,snrdb};
% rmse = RMSE_RRPRD{am,snrdb};
% elapsed = elapsed_time{am,snrdb};

for r = 1:100
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
AM2 = smooth(abs(cumsum(randn(size(t)))./Hz) + 1, 6000) ;
AM2 = 3*AM2 ./ max(AM2) - .7 ;
AM2 = AM2 .* exp(-((-8e4+1:2e4)'/5e4).^1.8) ;
ff = abs(cumsum(randn(1e5, 1))) ;   % standard brownian motion
% IF2 is shifted up from IF1 by exp(1)/2 +/- 1 ~ [0.35, 2.35]
IF2r = IF1r + ff./max(abs(ff)) + exp(1)/2 ;
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

%% signal samples
Lx = length(xn);
Tx = (0:length(xn)-1)./fs;

% % signal definition
% B = 400;
% phi1 = 100*Tx + B*Tx.^2;
% s1 = exp(2i*pi*phi1);
% sn = s1;
% 
% %% add complex noise
% noise = randn(1, Lx) + 1i*randn(1, Lx); % gaussian white noise
% 
% snr_in = -10;
% sn = add_noise(s1, noise, snr_in);

%% time frequency parameters
sigma_opt = 1/sqrt(600); % minimizing renyi entropy
[g, ~] = gauss_win(Lx, sigma_opt);

Nfft = length(xn)/2;
fr = fs/Nfft;

%% second order synchrosqueezing & ridge detection
fprintf("SST...");
[~, TFR] = sst2(xn, sigma_opt, Nfft);
% STFT = STFT(1:Nfft/2,:);

% P is the smoothness parameter for the weigthed spline approximation
% typicall value for synthetic signals at -10dB is 10^(-4)
tfrsq = TFR.SST2;
tfrsq(Nfft/2+1:end,:) = 0;
Fx = (0:Nfft-1)*fr;

P = 1 - 10^(-7);

fprintf("RD...\n");
% disp(elapsed(r))
try
    cALL = [];
    tic
    [Spl, ~] = RRP_RD(tfrsq, TFR.q_hat, TFR.omega1_hat, TFR.tau, P, 3);
    ii = 1;
    while ii <= 3
        cALL = [cALL ppval(Spl(ii).spline,(0:length(xn)-1)'./length(xn)) .* (fs/Nfft)];
        ii = ii + 1;
    end
    elapsed(r) = toc;
catch
    fprintf("ERR(3)!\n");
    try
        cALL = [];
        tic
        [Spl, ~] = RRP_RD(tfrsq, TFR.q_hat, TFR.omega1_hat, TFR.tau, P, 2);
        ii = 1;
        while ii <= 2
            cALL = [cALL ppval(Spl(ii).spline,(0:length(xn)-1)'./length(xn)) .* (fs/Nfft)];
            ii = ii + 1;
        end
        elapsed(r) = toc;
    catch
        fprintf("ERR(2)!\n");
        tic
        [Spl, ~] = RRP_RD(tfrsq, TFR.q_hat, TFR.omega1_hat, TFR.tau, P);
        cALL = ppval(Spl(1).spline,(0:length(xn)-1).'./length(xn)) .* (fs/Nfft);
        elapsed(r) = toc;
    end
end

cALL = cALL ./ 2;
[~,minIdx] = min(cALL(len/2,:));
c = cALL(:,minIdx);

curves(:,r) = c;
seg = 5*fs+1:35*fs;
rmse(r) = norm(c(seg)-IF1_fund(seg), 2)/norm(IF1_fund(seg), 2);
fprintf("%d, RMSE = %.2f, time = %.2f\n", r, rmse(r), elapsed(r));

curve_RRPRD{am,snrdb} = curves;
RMSE_RRPRD{am,snrdb} = rmse;
elapsed_time{am,snrdb} = elapsed;

save('./newARMAnoise/curve_RRPRD.mat', 'curve_RRPRD');
save('./newARMAnoise/RMSE_RRPRD.mat', 'RMSE_RRPRD');
save('./newARMAnoise/time_RRPRD.mat', 'elapsed_time');

end
end
end

%%
% fr = fs/Nfft;
% Fx = (0:Nfft-1)*fr;
% 
% figure;
% imageSQ(Tx, Fx, abs(tfrsq), 0.99); ylim([0 20])
% colormap(1-gray); colorbar;
% set(gca,'YDir','normal');
% xlabel("time");
% ylabel("frequency");
% title("short time Fourier transform");
% 
% figure;
% imageSQ(Tx, Fx, abs(tfrsq), 0.99); ylim([0 20]);
% colormap(1-gray); colorbar;
% set(gca,'YDir','normal');
% xlabel("time"); ylabel("frequency");
% hold on;
% plot(Tx, c, 'r');
% hold off;
% title("result of the ridge detection");
