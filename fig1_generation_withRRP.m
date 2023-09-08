% % % The overall flowchart
clear;
close all;
figPos = [100 50 550 550];
addpath('./TF_analysis/');
addpath(genpath('./MFEToolbox'));
addpath('./RRP-RD/toolbox/TF_Toolbox/');
addpath('./RRP-RD/RRP_alg/');

% for AFUND = [.0 1.0]
AFUND = [0.8 0.4 0.125 0.0625];
SNR = [Inf 20 10 5];
snrdb = 4;
am = 1;

% Generate signal sample
r = 1; %1:100
snr = SNR(snrdb);
disp([num2str(r),' / ',num2str(SNR(snrdb)),'db']);
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
AM1 = AM1 .* exp(-((t+10)/20).^2);

% generate FM for the 1st component
IF1r = smoothdata(abs(cumsum(randn(size(t)))), "gaussian", .4*length(t));
% IF1 is between pi/2+0.4+/-1; that is, around [1,3] Hz
IF1r = IF1r./max(abs(IF1r)) + 1.7;
% IF1 = smooth(IF1r, 8000);

% this is the smoothing step. Can use the shorter window to get more complicated IF
tmp = cumsum(exp(-((t-25)/2.5).^2));
tmp = 1.5*tmp ./ max(tmp);
IF1 = IF1r + tmp;

% phase
phi1 = [];
for l = 1:4
    u = randn(size(t))*sqrt(0.2);
    u = smoothdata(u, "gaussian", 0.1*length(t));
    phi1 = [phi1 l*(cumsum(IF1)./Hz + t.^1.9/(20*1.9)) + u];
end

% generate AM & FM for the 2nd component
AM2 = smooth(abs(cumsum(randn(size(t)))./Hz) + 1, 6000) ;
AM2 = 3*AM2 ./ max(AM2) - .5 ;
AM2 = AM2 .* exp(-((-1e4-1e5+1:0-1e4)'/4e4).^2);
% % IF2 is shifted up from IF1 by exp(1)/2 +/- 1 ~ [0.35, 2.35]
% ff = abs(cumsum(randn(1e5, 1))) ;
% IF2r = IF1r + ff./max(abs(ff)) + exp(1)/2 ;
% % this is the smoothing step. Can use the shorter window to get more complicated IF
% IF2 = smooth(IF2r, 5000) + t/20 ;
% phi2 = cumsum(IF2) ./ Hz + t.^2/(20*2) ;

AM = AM1 + AM2;
% figure; plot(AM1); hold on; plot(AM2); hold on;

% generate the 1st oscillatory component (like PPG)
aa = rand(1,3);
a4 = aa(3);
a3 = sum(aa(2:3)); 
a2 = sum(aa); a3 = 0.9*a2; a4 = 0.8*a3;
M = max([AFUND(am) a2]);
a1 = AFUND(am)./M;
a2 = a2./M;
a3 = a3./M;
a4 = a4./M;
s3 = (a1*cos(2*pi*phi1(:,1))) ...
    + (a2*cos(2*pi*phi1(:,2))) ...
    + (a3*cos(2*pi*phi1(:,3))) ...
    + (a4*cos(2*pi*phi1(:,4)));
x1 = AM .* s3;
hTrue = [AM.*(a1 * cos(2*pi*phi1(:,1))), ...
    AM.*(a2 * cos(2*pi*phi1(:,2))), ...
    AM.*(a3 * cos(2*pi*phi1(:,3))), ...
    AM.*(a4 * cos(2*pi*phi1(:,4)))];
hTrue = real(hTrue);
% figure; plot(hTrue(:,1));

x = x1;

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
y = random('T', 5, length(x)/2, 1);
% sigma = 10^(log10(std(x(1:length(x)/2))./std(y)) - (snr)/20);
noise(1:length(x)/2) = y;   % * sigma;

y = ARMA11(length(x)/2);
% sigma = 10^(log10(std(x(length(x)/2+1:end))./std(y)) - (snr)/20);
noise(length(x)/2+1:end) = y;   % * sigma;

sigma = 10^(log10(std(x)./std(noise)) - (snr)/20);
noise = noise * sigma;
fprintf("SNR = %.1f db\n", 20*log10(std(x)./std(noise)));
Y = x + noise;

seg = 26*Hz+1:34*Hz;
seg = 29*Hz+1:35*Hz;
figure;
plot(seg./Hz, real(Y(seg)),'LineWidth',1,'Color',[.7 .7 .7]);
hold on;
plot(seg./Hz, real(x(seg)),'LineWidth',1.8);
hold off;
ax = gca; ax.FontSize = 22; xlabel('time(sec)');
    
% downsample to 200Hz.
x = x(1:10:end) ;
x1 = x1(1:10:end);
Y = Y(1:10:end) ;
AM = AM(1:10:end) ;
IF1 = IF1(1:10:end) ;
% AM2 = AM2(1:10:end) ;
% IF2 = IF2(1:10:end) ;
t = t(1:10:end) ;
noise = noise(1:10:end) ;
orgHz = Hz;
Hz = Hz / 10 ;
hTrue = hTrue(1:10:end,:);

fs = 100;
% resample again to process faster
x = resample(x, fs, 200);
Y = resample(Y, fs, 200);
x1 = resample(x1, fs, 200);
hTrue = resample(hTrue, fs, 200);

IF1 = 2000*diff(phi1(:,1));
IF1 = resample(IF1, fs, 2000);
for k = 2:4
    tmp = 2000*diff(phi1(:,k));
    tmp = resample(tmp, fs, 2000);
    IF1 = [IF1 tmp];
end

%%
len = 2^(floor(log2(length(Y))));
xn = Y(1:len);
IF1_fund = IF1(1:len,1);

%% FMRD
sigma_opt = 1/sqrt(1000);
Nfft = len;
[~, ~, SST2, ~, ~, ~, ~, ~, tau2, tau3, phi22p, ~, ~] ...
    = sstn_test_new(xn, 0, sigma_opt, Nfft, 1:len);
% resample TFR
SST2 = SST2(1:2:end,:);
Nfft = Nfft / 2;
fr = fs/Nfft;
fr

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

tfrtic = (0:Nfft/2-1)*fr;
tfrtic = reshape(tfrtic, [], 1);

[~,minIdx] = min(cFM(1,:));
c = cFM(:,minIdx);
figure; set(gcf,'Position',figPos);
imageSQ((0:size(SST2,2)-1)./fs, tfrtic, abs(SST2)./max(max(abs(SST2))), 0.99);
axis xy; colormap(1-gray);
xlabel('time(sec)');
ylabel('frequency(Hz)');
set(gca, 'FontSize', 20);
ylim([0 25])
xlim([3 40])
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

figure; set(gcf,'Position',figPos);
imageSQ((0:size(tfrRRP,2)-1)./fs, FxRRP, abs(tfrRRP)./max(max(abs(tfrRRP))), 0.99);
ylim([0 25])
xlim([3 40])
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

%% basic parameters for STFT
basicTF.fr = 0.05; % frequency resolution
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)
% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here
advTF.Smo = 0; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection;
advTF.ths = 1E-9; % Global threshold of STFT
advTF.lpc = 0;

basicTF.hop = 1;
basicTF.win = fs*5+1;
basicTF.fs = fs; 
advTF.HighFreq = 50.0/basicTF.fs;   % highest frequency/sampling freq
advTF.LowFreq = 0.1/basicTF.fs;   % lowest frequency/sampling freq

clear P
% cepstrum
P.num_s = 1;
% (harmonic 的根數, 只做簡單探討的話設1就好)
P.num_c = 1;
% parameters for cepstral representation
cepR.g = 0.3; % for generalized cepstrum:
cepR.Tc = 0; % Global threshold of cepstrum

%% Get inverse cepstrum, STFT and SST
% r = 76;    
[tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(Y, advTF.LowFreq, advTF.HighFreq, ...
    basicTF.fr/fs, basicTF.hop, basicTF.win, 1, 5, 1, 1, 0);
[h, ~, ~] = hermf(basicTF.win,1,5);
h0 = h(floor(size(h,2)/2)+1);

% % TFR plot
% TFR = dsSST; TFRtic = dstfrtic;
% figure; set(gcf,'Position',figPos);
% imageSQ((0:size(TFR,2)-1)./fs, TFRtic*fs, abs(TFR), 0.99);
% axis xy; colormap(1-gray); %colorbar
% xlabel('time(sec)'); ylabel('frequency(Hz)');
% ax = gca; ax.FontSize = 20;
% ylim([0 17])
% xlim([3 40])
% 
% figure; set(gcf,'Position',figPos);
% imageSQ((0:size(TFR,2)-1)./fs, TFRtic*fs, abs(TFR), 0.99);
% axis xy; colormap(1-gray); %colorbar
% xlabel('time(sec)','FontSize',20);
% % ylabel('frequency(Hz)','FontSize',20);
% % title(['D = ',num2str(AF)]);
% ax = gca; ax.FontSize = 20; ylim([0 12]);
% hold on;
% % plot((0:length(IF)-1)./fs, IF, '--b', 'LineWidth', 4);
% plot((0:length(c)-1)./fs, TFRtic(c)*fs, 'r', 'LineWidth', 4);
% hold off;
% ylabel('frequency(Hz)')
% ylim([0 17])
% xlim([3 40])

figure; set(gcf,'Position',figPos);
imageSQ((0:size(tfr,2)-1)./fs, tfrtic*fs, abs(tfr)./max(max(abs(tfr))), 0.99);
axis xy; colormap(1-gray); %colorbar
xlabel('time(sec)','FontSize',20);
ylabel('frequency(Hz)','FontSize',20);
ax = gca; ax.FontSize = 20;
ylim([0 25])
xlim([3 40])

figure; set(gcf,'Position',figPos);
imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(tfrsq)./max(max(abs(tfrsq))), 0.99);
axis xy; colormap(1-gray);
xlabel('time(sec)','FontSize',20);
ylabel('frequency(Hz)','FontSize',20);
ax = gca; ax.FontSize = 20;
ylim([0 25])
xlim([3 40])

%% SingleCurveExt; peeling singleCurveExt
cSH = [];
TFR = abs(tfrsq);
for k = 1:4
    disp(['curve ',num2str(k)]);
    idx = find(tfrsqtic*fs > 0.5 & tfrsqtic*fs < 20);
    [c] = CurveExt(TFR(idx,:).', 1);
    c = idx(1) + c - 1;
    cSH = [cSH reshape(c, [], 1)];
    for t = 1:length(c)
        TFR(c(t)-round(0.3/basicTF.fr) : c(t)+round(0.3/basicTF.fr), t) = 0;
    end
end

[~,minIdx] = min(cSH(1,:));
c = cSH(:,minIdx);

figure; set(gcf,'Position',figPos);
imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
axis xy; colormap(1-gray); %colorbar
xlabel('time(sec)','FontSize',20);
set(gca, 'FontSize', 20);
ax.YTickLabel = [];
hold on;
plot((0:size(cSH,1)-1)./fs, tfrsqtic(c)*fs, 'r', 'LineWidth', 3.0);
hold on; plot((0:length(IF1_fund)-1)./fs, IF1_fund, '--b', 'LineWidth', 4.0);
% for k = 1:4
%     IF = 2000*diff(k*phi1);
%     IF = resample(IF, basicTF.fs, 2000);
%     % hold on; plot((0:size(IF,1)-1)./fs, IF, 'b', 'LineWidth', 1.0);
%     hold on; plot((0:size(cSH,1)-1)./fs, tfrsqtic(cSH(:,k))*fs, 'r', 'LineWidth', 3.0);
% end
% hold on; plot((0:length(IF1)-1)./fs, IF1, '--b', 'LineWidth', 6); hold off;
ylim([0 25])
xlim([3 40])

%% MHRD
cMH = [];
cEND = [1; 1; 1];
lambda = [1.0 0.8 0.7];
mu = [.0 .0];
num_multi = [3 3];
band_fund = [0.5 8.0];
band_multi = [0.8 0.5];
ver = 3;
bw = [-0.8 0.8];

tt = [3*fs:2*fs:size(tfrsq,2)];
tt = [tt size(tfrsq,2)];

% MHRD(2)
cMH = MultiCurveExt(tfrsq, tfrsqtic*fs, lambda(1:2), mu(1), num_multi(1), ...
    band_fund, band_multi(1), 2, tt, bw);
cMH2 = CurveMultiExt_withFund(abs(tfrsq(:,3*fs+1:end)).', tfrsqtic, cMH(:,1), 0.7, ...
    0, 2, round(1.0/basicTF.fr));
cMH4 = CurveMultiExt_withFund(abs(tfrsq(:,3*fs+1:end)).', tfrsqtic, cMH(:,1), 0.7, ...
    0, 4, round(1.0/basicTF.fr));

cMH = [cMH(:,1) cMH2 cMH(:,2) cMH4];

% % MHRD(3)
% cMH = MultiCurveExt(tfrsq, tfrsqtic*fs, lambda, mu, num_multi, ...
%     band_fund, band_multi, ver, tt, bw);
% [cMH(:,4)] = CurveMultiExt_withFund(abs(tfrsq(:,3*fs+1:end)).', tfrsqtic, cMH(:,1), 0.6, ...
%     0, 4, round(1/basicTF.fr));

% padding ones
cMH = [ones(3*fs,size(cMH,2))*diag(cMH(1,:)); cMH];

figure; set(gcf,'Position',figPos);
imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
axis xy; colormap(1-gray); %colorbar
xlabel('time(sec)','FontSize',20);
% ylabel('frequency(Hz)','FontSize',20);
ax = gca; ax.FontSize = 20;
for k = 1:size(cMH,2)
    hold on; plot((0:size(cMH,1)-1)./fs, tfrsqtic(cMH(:,k)).*fs, '-r', 'LineWidth', 3);
    if k==1
    hold on; plot((0:length(IF1_fund)-1)./fs, IF1_fund, '--b', 'LineWidth', 4.0);
    end
end
hold off;
ylim([0 25])
xlim([3 40])

%% Recon
hw = [];
bw = [0.4 0.6 0.8 1];
for k = 1:4
    IF = 2000*diff(k*phi1);
    IF = resample(IF, basicTF.fs, 2000);
    c_true = zeros(size(IF));
    for t = 1:length(c_true)
        [~,c_true(t)] = min(abs(tfrsqtic*fs-IF(t)));
    end
    tmp = Recon_sqSTFT_v2(tfrsq, tfrsqtic, fs, cMH(:,k), bw(k), h0);
    hw = [hw reshape(tmp,[],1)];
end

seg = 29*fs+1:35*fs;
% seg = 28*fs+1:32.5*fs;

figure;
for k = 1:4
    subplot(4,1,4-k+1);
    plot(seg./fs,real(hw(seg,k)),'-r','LineWidth',2); hold on;
    plot(seg./fs,hTrue(seg,k),'-k','LineWidth',.8); hold off;
    ax = gca; ax.FontSize = 22;
    if k == 1; xlabel('time(sec)'); else; ax.XTickLabel = []; end
    xlim([29 35])
    % if k==2; ylim([-2 2]); else; ylim([-1 1]); end
end
hwALL = real(sum(hw.')).';

figure;
plot(seg./fs, hwALL(seg), 'r', 'LineWidth', 2); hold on;
plot(seg./fs, real(x(seg)), 'k', 'LineWidth', 1); hold off;
ax = gca; ax.FontSize = 22; xlabel('time(sec)');
xlim([29 35])
