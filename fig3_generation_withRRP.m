% % % "RD on two IMT components demo" figure
clear;
close all;
addpath('./TF_analysis/');
addpath(genpath('./MFEToolbox/'));
addpath('./RRP-RD/toolbox/TF_Toolbox/');
addpath('./RRP-RD/RRP_alg/');

AFUND1 = 0.75;
AFUND2 = 2.0;
SNR = [Inf 20 10 5 0];

snrdb = 1;
% Generate signal sample
r = 76; %1:100
disp([num2str(r),' / ',num2str(SNR(snrdb)),'db']);
    % change this random seed if you want to see other results
rand('state', r);
randn('state', r);

%% generate simulated signals 
% start with a high sampling rate so that the phase can be more
% accurately generated (will downsample it later)
Hz = 2000 ;
t = (1:1e5)' / Hz ;

% generate AM & FM for the 1st component
AM1 = smooth(abs(cumsum(randn(size(t)))./Hz) + 1, 2500) ;
AM1 = 3*AM1 ./ max(AM1) - .5 ; 
AM1 = AM1 .* exp(-((-2e4+1:8e4)'/5e4).^2) ;
ff = abs(cumsum(randn(size(t)))) ;
% IF1 is between pi/2+0.4+/-1; that is, around [1,3] Hz
IF1r = ff./max(abs(ff)) + (pi/2+0.4-1) ;
% this is the smoothing step. Can use the shorter window to get more
% complicated IF
tmp = cumsum(exp(-((t-25)/2.5).^2)) ; tmp = 1.5*tmp ./ max(tmp) ;
IF1 = smooth(IF1r, 8000) + tmp ;
phi1 = cumsum(IF1) ./ Hz + t.^1.9/(20*1.9) ;

% generate AM & FM for the 2nd component
AM2 = smooth(abs(cumsum(randn(size(t)))./Hz) + 1, 6000) ;
AM2 = AM2 ./ max(AM2) - .7 ;
AM2 = 3 * AM2 .* exp(-((-8e4+1:2e4)'/5e5).^1.8) ;
ff = abs(cumsum(randn(1e5, 1))) ;   % standard brownian motion
% IF2 is shifted up from IF1 by exp(1)/2 +/- 1 ~ [0.35, 2.35]
IF2r = IF1r + ff./max(abs(ff)) + exp(1)/2 ;
% this is the smoothing step. Can use the shorter window to get more
% complicated IF
IF2 = smooth(IF2r, 5000) + t/20 ;   % smooth(standard brownian) + e/2 + t/20
phi2 = cumsum(IF2) ./ Hz + t.^2/(20*2) ;

% generate the 1st oscillatory component (like PPG)
aa = rand(1,3); a4 = aa(1); a3 = sum(aa(1:2)); a2 = sum(aa);
M = max([AFUND1 a2]);
a1 = AFUND1./M;
a2 = a2./M;
a3 = a3./M;
a4 = a4./M;
ph = rand(1,4)*2*pi-pi;
s3 = a1.*cos(2*pi*phi1+ph(1)) ...
    + a2.*cos(2*pi*2*phi1+ph(2)) ...
    + a3.*cos(2*pi*3*phi1+ph(3)) ...
    + a4.*cos(2*pi*4*phi1+ph(4));
x1 = AM1 .* s3;
x1_mode1 = AM1 .* (a1 * cos(2*pi*phi1+ph(1)));

% generate the 2nd oscillatory component (with the very general shape) you can change the shape to contain only 5 or 8 harmonics
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

% generate the 2nd oscillatory component (like PPG)
bb = rand(1,4); b4 = bb(4); b3 = sum(bb(3:4)); b2 = sum(bb(2:4)); b1 = sum(bb);
M = max([b1 b2 b3 b4]);
b1 = AFUND2*b1./M ;
b2 = AFUND2*b2./M ;
b3 = AFUND2*b3./M;
ph = rand(1,4)*2*pi-pi;
s3 = b1.*cos(2*pi*phi2+ph(1)) ...
    + b2.*cos(2*pi*2*phi2+ph(2)) ...
    + b3.*cos(2*pi*3*phi2+ph(3)) ...
    + 0.6*b3.*cos(2*pi*4*phi2+ph(4)); 
x2 = AM2 .* s3;
x2_mode1 = AM2 .* (b1 * cos(2*pi*2*phi1+ph(1)));

x = x1 + x2 ;

noise = random('T',5,length(x2),1) ;
sigma = 10^(log10(std(x)./std(noise)) - SNR(snrdb)/20);
noise = sigma * noise ; 
% var(noise)
% snrdb = 20 * log10(std(x)./std(noise)) ;
fprintf(['SNR = ',num2str(20*log10(std(x)./std(noise))),'\n']);

Y = x + noise;

% downsample to 200Hz.
x = x(1:10:end) ;
x1 = x1(1:10:end);
x2 = x2(1:10:end);
Y = Y(1:10:end) ;
AM1 = AM1(1:10:end) ;
IF1 = IF1(1:10:end) ;
AM2 = AM2(1:10:end) ;
IF2 = IF2(1:10:end) ;
t = t(1:10:end) ;
noise = noise(1:10:end) ;
orgHz = Hz;
Hz = Hz / 10;
x1_mode1 = x1_mode1(1:10:end);
x2_mode1 = x2_mode1(1:10:end);

fs = 100;
x = resample(x, fs, 200);  % resample again to process faster
Y = resample(Y, fs, 200);
x1 = resample(x1, fs, 200);
x2 = resample(x2, fs, 200);
x1_mode1 = resample(x1_mode1, fs, 200);
x2_mode1 = resample(x2_mode1, fs, 200);

IF = 2000*diff(phi2); IF = resample(IF, fs, 2000);

IF1 = 2000*diff(phi1); IF1 = resample(IF1, fs, 2000);
for k = 2:4
    tmp = 2000*diff(k*phi1); tmp = resample(tmp, fs, 2000);
    IF1 = [IF1 tmp];
end

len = 2^(floor(log2(length(Y))));
xn = Y(1:len);
IF1_fund = IF1(1:len,1);

%% FMRD
sigma_opt = 1/sqrt(1000);
Nfft = len;
[STFT, ~, SST2, ~, ~, ~, ~, ~, tau2, tau3, phi22p, ~, ~] ...
    = sstn_test_new(xn, 0, sigma_opt, Nfft, 1:len);
% resample TFR
SST2 = SST2(1:2:end,:);
Nfft = Nfft / 2;
fr = fs/Nfft;
fr
tfrtic = (0:Nfft/2-1)*fr;
tfrtic = reshape(tfrtic, [], 1);

%%
cFM = [];
TFR = SST2;
% peeling for RD
fprintf("RD: ");
bw = round(0.2/fr);
tic
for k = 1:4
    fprintf("%d ... ", k);
    q = ceil(abs(Nfft*real(phi22p)/len^2));
    c = exridge_CR_MB(TFR, 0, 0, q, 3);
    cFM = [cFM reshape(c,[],1)];
    for l = 1:length(c)
        % masking
        TFR(max(c(l)-bw,1):min(c(l)+bw,size(TFR,1)) , l) = 0;
    end
end
fprintf("\n");

[~,minIdx] = min(cFM(1,:));
c = cFM(:,minIdx);
figure;
% set(gcf,'Position',figPos);
imageSQ((0:size(SST2,2)-1)./fs, tfrtic, abs(SST2)./max(max(abs(SST2))), 0.99);
axis xy; colormap(1-gray);
xlabel('time(sec)');
ylabel('frequency(Hz)');
set(gca, 'FontSize', 20);
ylim([0 20])
xlim([10 40])
% hold on;
% plot((0:size(cFM,1)-1)./fs, tfrtic(c), 'r', 'LineWidth', 1);
for i = 1:4
    hold on;
    plot((1:length(cFM(:,i)))./fs, tfrtic(cFM(:,i)), 'r', 'LineWidth', 3);
end
hold off;

%% RRP-RD
sigma_opt = 1/sqrt(1000); % minimizing renyi entropy
[g, ~] = gauss_win(len, sigma_opt);

% second order synchrosqueezing & ridge detection
Nfft = length(xn)/2;
fr = fs/Nfft
[~,TFR] = sst2(xn, sigma_opt, Nfft);
% tfrsq for RRP-RD
tfrRRP = TFR.SST2;
tfrRRP(Nfft/2+1:end,:) = 0;
FxRRP = (0:Nfft-1)*fr;

P = 1 - 10^(-7);

fprintf("RD...\n");
% disp(elapsed(r))
try
    cRRP = [];
    tic
    [Spl, ~] = RRP_RD(tfrRRP, TFR.q_hat, TFR.omega1_hat, TFR.tau, P, 4);
    ii = 1;
    while ii <= 4
        cRRP = [cRRP ppval(Spl(ii).spline,(0:length(xn)-1)'./length(xn)) .* (fs/Nfft)];
        ii = ii + 1;
    end
    elapsed(r) = toc;
catch
    fprintf("ERR(4)!\n");
    try
        cRRP = [];
        tic
        [Spl, ~] = RRP_RD(tfrRRP, TFR.q_hat, TFR.omega1_hat, TFR.tau, P, 3);
        ii = 1;
        while ii <= 3
            cRRP = [cRRP ppval(Spl(ii).spline,(0:length(xn)-1)'./length(xn)) .* (fs/Nfft)];
            ii = ii + 1;
        end
        elapsed(r) = toc;
    catch
        fprintf("ERR(3)!\n");
        try
            cRRP = [];
            tic
            [Spl, ~] = RRP_RD(tfrRRP, TFR.q_hat, TFR.omega1_hat, TFR.tau, P, 2);
            ii = 1;
            while ii <= 2
                cRRP = [cRRP ppval(Spl(ii).spline,(0:length(xn)-1)'./length(xn)) .* (fs/Nfft)];
                ii = ii + 1;
            end
            elapsed(r) = toc;
        catch
            fprintf("ERR(2)!\n");
            tic
            [Spl, ~] = RRP_RD(tfrRRP, TFR.q_hat, TFR.omega1_hat, TFR.tau, P);
            cRRP = ppval(Spl(1).spline,(0:length(xn)-1).'./length(xn)) .* (fs/Nfft);
            elapsed(r) = toc;
        end
    end
end

cRRP = cRRP ./ 2;
[~,minIdx] = min(cRRP(len/2,:));
c = cRRP(:,minIdx);

figure;
imageSQ((0:size(tfrRRP,2)-1)./fs, FxRRP, abs(tfrRRP)./max(max(abs(tfrRRP))), 0.99);
ylim([0 20])
xlim([10 40])
colormap(1-gray);
set(gca, 'FontSize', 20);
xlabel("time");
ylabel("frequency");
for i = 1:size(cRRP,2)
    hold on;
    plot((0:size(cRRP,1)-1)./fs, cRRP(:,i), '-r', 'LineWidth', 3.0);
end
% hold on;
% plot((0:size(cRRP,1)-1)./fs, cRRP(:,1), '-r', 'LineWidth', 3.0);
hold on;
plot((0:length(IF1_fund)-1)./fs, IF1_fund, '--b', 'LineWidth', 4.0);
hold off;

%%
% %% basic parameters for STFT
% basicTF.fr = fr; % frequency resolution
% basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)
% % advanced parameters for STFT
% advTF.num_tap = 1; % Number of tap in ConceFT
% advTF.win_type = 'Gauss'; % Only 2-tap basis here
% advTF.Smo = 0; % alpha smoothing; 1 = no smoothing
% advTF.Rej = 0; % The bandwidth of window rejection;
% advTF.ths = 1E-9; % Global threshold of STFT
% advTF.lpc = 0;
% 
% basicTF.hop = 1;
% basicTF.win = fs*5+1;
% basicTF.fs = fs; 
% advTF.HighFreq = 50.0/basicTF.fs;   % highest frequency/sampling freq
% advTF.LowFreq = 0.01/basicTF.fs;   % lowest frequency/sampling freq
% 
% clear P
% % cepstrum
% P.num_s = 1;
% % (harmonic 的根數, 只做簡單探討的話設1就好)
% P.num_c = 1;
% % parameters for cepstral representation
% cepR.g = 0.3; % for generalized cepstrum:
% cepR.Tc = 0; % Global threshold of cepstrum

segTFR = [10 40];
segsig = [20 30];

figure; plot((1:length(x1))./fs, x1, 'k', 'LineWidth', 1.5); xlim(segsig);
set(gca, 'FontSize', 20);
figure; plot((1:length(x2))./fs, real(x2), 'k', 'LineWidth', 1.5); xlim(segsig);
set(gca, 'FontSize', 20); xlabel('time(sec)');

%% Get inverse cepstrum, STFT and SST
% % SST
% r = 76;
% [h, ~, ~] = hermf(basicTF.win,1,5);
% h0 = h(floor(size(h,2)/2)+1);
% [tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(Y, advTF.LowFreq, advTF.HighFreq, ...
%     basicTF.fr/fs, basicTF.hop, basicTF.win, 1, 5, 1, 1, 0);

tfrsq = TFR.SST2;
tfr = STFT;
% tfrsqtic = tfrtic;
tfrsqtic = reshape(FxRRP, [], 1);

figure;
subplot(1,3,2);imageSQ((1:size(tfrsq,2))./fs, tfrsqtic, abs(tfrsq), 0.99);
axis xy; colormap(1-gray); %colorbar
xlabel('time(sec)');
% ylabel('frequency(Hz)');
set(gca, 'FontSize', 20);
ylim([0 20]);
xlim(segTFR);

subplot(1,3,1);
imageSQ((1:size(tfr,2))./fs, tfrtic, abs(tfr), 0.99);
axis xy; colormap(1-gray); %colorbar
xlabel('time(sec)');
ylabel('frequency(Hz)');
set(gca, 'FontSize', 20);
ylim([0 20]);
xlim(segTFR);

% tfrtic = (0:Nfft/2-1)*fr;
% tfrtic = reshape(tfrtic, [], 1);
[~,minIdx] = min(cRRP(1,:));
c = cFM(:,minIdx);
subplot(1,3,3);
imageSQ((1:size(tfrsq,2))./fs, tfrsqtic, abs(tfrsq), 0.99);
for k = 1:3
    hold on;
    plot((1:size(cALL,1))./fs, tfrsqtic(cALL(:,k)), '-r', 'LineWidth', 3);
end
hold on;
plot((1:length(c))./fs, tfrtic(c), ':b', 'LineWidth', 3);
% for k = 1:3
%     hold on;
%     plot((1:length(c))./fs, cRRP(:,k), '-b', 'LineWidth', 1.2);
% end
axis xy; colormap(1-gray); %colorbar
xlabel('time(sec)');
% ylabel('frequency(Hz)');
set(gca, 'FontSize', 20);
ylim([0 20]);
xlim(segTFR);

%% 3curves-MultiCurveExt
lambda = [1.0, 0.8, 0.6];
mu = [0.3 0];
multi = [2, 3];
band = [0.3, 8.0];
band_multi = [.2, .4];
cALL = [];
tt = 0:fs/2:size(tfrsq,2);

cALL = MultiCurveExt(tfrsq, tfrsqtic, lambda, mu, multi, ...
    band, band_multi, 3, tt, [-.3 .5]);

figure;
imageSQ((1:size(tfrsq,2))./fs, tfrsqtic, abs(tfrsq), 0.99);
axis xy; colormap(1-gray); colorbar
xlabel('time(sec)');
% ylabel('frequency(Hz)','FontSize',20);
ax = gca; ax.FontSize = 20; ylim([0 20]);
for k = 1:3
    hold on;
    plot((1:size(cALL,1))./fs, tfrsqtic(cALL(:,k)), '--r', 'LineWidth', 3);
end
hold off;
xlim(segTFR);

% %% comb filter
% comb = abs(tfrsq);
% for t = 1:size(comb,2)
% for f = 1:size(comb,1)
% for k = 2:4
%     RR = max(1, k*f-10) : min(size(comb,1), k*f+10);
%     add = max(abs(tfrsq(RR,t)));
%     if isempty(RR)
%         add = 0;
%     end
%     comb(f,t) = comb(f,t) + add;
% end
% end
% end

%%
figure;
imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(comb), 0.99);
colormap(1-gray); colorbar;
set(gca,'YDir','normal');
xlabel("time"); ylabel("frequency"); ylim([0 15])
title("result of the ridge detection");

idx = find(tfrsqtic*fs > 0.5 & tfrsqtic*fs < 10);
[c] = CurveExt(abs(comb(idx,:)).', .0);
c = c + idx(1) - 1;
hold on;
plot((0:length(c)-1)./fs, tfrsqtic(c)*fs, 'r');
hold off;
