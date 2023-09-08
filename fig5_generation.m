% % % figure of SAMD example
clear;
close all;
addpath('./TF_analysis/');
addpath(genpath('./MFEToolbox/'));
addpath('./RRP-RD/toolbox/TF_Toolbox/');
addpath('./RRP-RD/RRP_alg/');

AFUND1 = 1.5; AFUND2 = 1.0; SNR = [Inf 20 10 5 0];

snrdb = 4;
snr = SNR(snrdb);

D = AFUND1;
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
t = (1:1e5).' / Hz ;

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
for l = 1:4
    u = randn(size(t))*sqrt(0.1);
    u = smoothdata(u, "gaussian", 0.1*length(t));
    phi1 = [phi1 l*(cumsum(IF1)./Hz + t.^1.9/(20*1.9)) + u];
end

% generate the 1st oscillatory component (like PPG)
aa = rand(1,2); %a4 = aa(1);
a3 = aa(1); a4 = a3 * 0.8;
a2 = sum(aa);
M = max([D a2 a3]);
a1 = D./M; a2 = a2./M; a3 = a3./M; a4 = a4./M;
s3 = a1 * cos(2*pi*phi1(:,1)) ...
    + a2 * cos(2*pi*phi1(:,2)) ...
    + a3 * cos(2*pi*phi1(:,3)) ...
    + a4 * cos(2*pi*phi1(:,4));
x1 = AM1 .* s3;
x1_mode1 = AM1 .* (a1 * cos(2*pi*phi1(:,1)));

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
% generate AM & FM for the 2nd component

AM2 = smooth(abs(cumsum(randn(size(t)))./Hz) + 1, 6000) ;
AM2 = 3*AM2 ./ max(AM2) - .7 ;
AM2 = AM2 .* exp(-((-8e4+1:2e4)'/5e5).^1.8) ;
ff = abs(cumsum(randn(1e5, 1))) ;   % standard brownian motion
% IF2 is shifted up from IF1 by exp(1)/2 +/- 1 ~ [0.35, 2.35]
IF2r = IF1r + ff./max(abs(ff)) + exp(1)/2 ;
% this is the smoothing step. Can use the shorter window to get more
% complicated IF
IF2 = smooth(IF2r, 5000) + t/20 ;   % smooth(standard brownian) + e/2 + t/20
phi2 = [];
for l = 1:4
    u = randn(size(t))*sqrt(0.1);
    u = smoothdata(u, "gaussian", 5000);
    phi2 = [phi2 l*(cumsum(IF2)./Hz + t.^2/(2*20)) + u];
end

% generate the 2nd oscillatory component (like PPG)
bb = rand(1,3);
b4 = bb(1);
b3 = sum(bb(1:2));
b2 = sum(bb);
M = max([AFUND2 b2]);
b1 = AFUND2./M;
b2 = b2./M;
b3 = b3./M;
b4 = b4./M;
s3 = b1.*cos(2*pi*phi2(:,1)) ...
    + b2.*cos(2*pi*phi2(:,2)) ...
    + b3.*cos(2*pi*phi2(:,3)) ...
    + b4.*cos(2*pi*phi2(:,4)); 
x2 = AM2 .* s3;
x2_mode1 = AM2 .* (b1.*cos(2*pi*phi2(:,1)));

x = x1 + x2 ;

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


seg = 25*Hz+1:30*Hz;
figure; plot(seg./Hz, real(Y(seg)),'LineWidth',1,'Color',[.7 .7 .7]); hold on; plot(seg./Hz, real(x(seg)),'LineWidth',1.5); hold off;
ax = gca; ax.FontSize = 24; xlabel('time(sec)')

% downsample to 200Hz.
res = Hz / 200;
x = x(1:res:end) ;
x1 = x1(1:res:end);
x2 = x2(1:res:end);
Y = Y(1:res:end) ;
AM1 = AM1(1:res:end) ;
IF1 = IF1(1:res:end) ;
AM2 = AM2(1:res:end) ;
IF2 = IF2(1:res:end) ;
t = t(1:res:end) ;
noise = noise(1:res:end) ;
orgHz = Hz;
Hz = Hz / res;
x1_mode1 = x1_mode1(1:res:end);
x2_mode1 = x2_mode1(1:res:end);

% x = resample(x, 100, 200);  % resample again to process faster
% Y = resample(Y, 100, 200);
% noise = resample(noise, 100, 200);
% x1 = resample(x1, 100, 200);
% x2 = resample(x2, 100, 200);
% hTrue1 = resample(hTrue1, 100, 200);

fs = 200;
IFtrue2 = 2000*diff(phi2(:,1)); 
IFtrue2 = resample(IFtrue2, fs, 2000);

IFtrue1 = 2000*diff(phi1(:,1));
IFtrue1 = resample(IFtrue1, fs, 2000);
for k = 2:4
    tmp = 2000*diff(phi1(:,k));
        tmp = resample(tmp, fs, 2000);
    IFtrue1 = [IFtrue1 tmp];
end

seg = 25*fs+1:30*fs;

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

% cepstrum
P.num_s = 1; % (harmonic 的根數, 只做簡單探討的話設1就好)
P.num_c = 1;
% parameters for cepstral representation
cepR.g = 0.3; % for generalized cepstrum:
cepR.Tc = 0; % Global threshold of cepstrum

%% Get inverse cepstrum, STFT and SST
r = 76;
[h, ~, ~] = hermf(basicTF.win,1,5);
h0 = h(floor(size(h,2)/2)+1);
[tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(Y, advTF.LowFreq, advTF.HighFreq, ...
    basicTF.fr/fs, basicTF.hop, basicTF.win, 1, 5, 1, 1, 0);

% figure; set(gcf,'Position',[100 120 1000 450]);
% imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
% axis xy; colormap(1-gray); colorbar
% xlabel('time(sec)');
% ylabel('frequency(Hz)');
% ax = gca; ax.FontSize = 20; ylim([0 25]);

%%
figure;
set(gcf,'Position',[100 120 1000 450]);
subplot(1,3,1); imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
axis xy; colormap(1-gray); colorbar
xlabel('time(sec)');
ylabel('frequency(Hz)');
ax = gca; ax.FontSize = 20; ylim([0 25]); xlim([20 30])

subplot(1,3,2);
imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
hold on;
for k = 1:4
    plot((0:size(IFtrue1,1)-1)./fs, IFtrue1(:,k), '--b', 'LineWidth', 3);
end
axis xy; colormap(1-gray); colorbar
xlabel('time(sec)');
ylabel('frequency(Hz)');
ax = gca; ax.FontSize = 20; ylim([0 25]); xlim([20 30])

%% 3curves
cALL = [];
cEND = [1; 1; 1];
lambda = [1.0, 0.8, 0.6];
mu = [0.3 0.2];
multi = [2, 3];
band = [1.0, 6.0];
band_multi = [.2, .4];
bw = [-0.3 0.5];
tt = 0:fs:size(tfrsq,2);

cALL = MultiCurveExt(tfrsq, tfrsqtic*fs, lambda, mu, multi, ...
    band, band_multi, 3, tt, bw);

% figure;
% set(gcf,'Position',[100 50 1000 700]);
% imageSQ((1:size(tfrsq,2))./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
% axis xy; colormap(1-gray); colorbar
% xlabel('time(sec)','FontSize',20);
% ylabel('frequency(Hz)','FontSize',20);
% ax = gca; ax.FontSize = 20; ylim([0 20]);
% for k = 1:3
%     hold on; plot((1:size(cALL,1))./fs, tfrsqtic(cALL(:,k)).*fs, 'r', 'LineWidth',1);
% end
% hold off;

%%
idx = find(tfrsqtic*fs > 0.1 & tfrsqtic*fs < 6.0);
[c] = CurveExt(abs(tfrsq(idx,:)'), 1.0);
c = idx(1) - 1 + c;

% figure;
% set(gcf,'Position',[100 50 1000 700]);
% imageSQ((1:size(tfrsq,2))./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
% axis xy; colormap(1-gray); colorbar
% xlabel('time(sec)','FontSize',20);
% ylabel('frequency(Hz)','FontSize',20);
% ax = gca; ax.FontSize = 20; ylim([0 20]);
% hold on;
% plot((1:size(c,1))./fs, tfrsqtic(c).*fs, 'r', 'LineWidth',1);
% hold off;

%% SAMD
% RECONSTRUCTION of the fundamental
f1_fund = Recon_sqSTFT_v2(tfrsq, tfrsqtic, fs, cALL(:,1), 0.3, h0).';
phi1_est = unwrap(angle(f1_fund));

% Fundamental of 2nd extraction and reconstruction???
% idx = find(tfrsqtic*fs > 0.3 & tfrsqtic*fs < 10.0);
% [c] = CurveExt(abs(tfrsq(idx,:)'), 1.0);
% c = idx(1) - 1 + c;
% seg = 23*fs+1:27*fs;
% disp(['IF1_error = ', num2str(norm(IF1(seg)-tfrsqtic(cALL(seg,1))*fs)/norm(IF1(seg)))]);
% disp(['fund1_error = ', num2str(norm(hTrue1(seg)-real(f1_fund(seg)))/norm(hTrue1(seg)))]);
% hold on;
% plot((0:length(c)-1)./fs, tfrsqtic(c)*fs, 'r');
% hold off;

AM1_est = abs(f1_fund).';
ii = 1;

Y1 = Y;

I = 3;
while ii <= I
[~, f1est, ~] = SAMD(real(Y1).', AM1_est, phi1_est  .', 4, 2);

% The 2nd IMT component
Y2 = real(Y) - f1est.';
Y2 = Y2 - mean(Y2);

% SST
[~, ~, tfrsq2, ~, tfrsqtic2] = ConceFT_sqSTFT_C(Y2, advTF.LowFreq, advTF.HighFreq, ...
    basicTF.fr/fs, basicTF.hop, basicTF.win, 1, 5, 1, 1, 0);
[h, ~, ~] = hermf(basicTF.win,1,5);
h0 = h(floor(size(h,2)/2)+1);
% 2nd IMT detection
tt = [0 3*fs:fs:size(tfrsq2,2)];
c2 = MultiCurveExt(tfrsq2, tfrsqtic2*fs, [1.0 0.8], 0, 2, ...
    [0.3 10.0], 0.3, 2, tt, [-1.0 1.0]);
f2_fund = Recon_sqSTFT_v2(tfrsq2, tfrsqtic2, fs, c2(:,1), 0.2, h0).';
phi2_est = unwrap(angle(f2_fund));
AM2_est = abs(f2_fund).';

% SAMD
[~, f2est, ~] = SAMD(Y2.', AM2_est, phi2_est.', 4, 2);

if ii==I
    subplot(1,3,3);
    imageSQ((0:size(tfrsq2,2)-1)./fs, tfrsqtic2*fs, abs(tfrsq2), 0.99);
    axis xy; colormap(1-gray); colorbar
    xlabel('time(sec)','FontSize',20);
    % ylabel('frequency(Hz)','FontSize',20);
    ax = gca; ax.FontSize = 20; ylim([0 25]); xlim([20 30])
end

Y1 = Y - f2est.';
Y1 = Y1 - mean(Y1);
[~, ~, tfrsq2, ~, tfrsqtic2] = ConceFT_sqSTFT_C(Y1, advTF.LowFreq, advTF.HighFreq, ...
    basicTF.fr/fs, basicTF.hop, basicTF.win, 1, 5, 1, 1, 0);
[h, ~, ~] = hermf(basicTF.win,1,5);
h0 = h(floor(size(h,2)/2)+1);

% Usually we don't need to estimate the fundamental IF for the 1st IMT again
% idx = find(tfrsqtic2*fs > 0.3 & tfrsqtic2*fs < 10.0);
% [c1] = CurveExt(abs(tfrsq2(idx,:)'), 1.0);
% c1 = idx(1) - 1 + c1;

f1_fund = Recon_sqSTFT_v2(tfrsq2, tfrsqtic2, fs, cALL(:,1), 0.4, h0).';
phi1_est = unwrap(angle(f1_fund));
AM1_est = abs(f1_fund).';

ii = ii + 1;

end

%%
% seg = 30*fs+1:35*fs;
fest = f1est + f2est;
figure;
plot(seg./fs, real(f1est(seg)),'-r','LineWidth',2.5);
hold on;
plot(seg./fs, real(x1(seg)),'-k','LineWidth',1.2);
hold off;
ax = gca; ax.FontSize = 24; ax.XTickLabel = [];

figure;
plot(seg./fs, real(f2est(seg)),'-r','LineWidth',2.5);
hold on;
plot(seg./fs, real(x2(seg)),'-k','LineWidth',1.2);
hold off;
ax = gca; ax.FontSize = 24; ax.XTickLabel = [];

figure;
plot(seg./fs, real(fest(seg)),'-r','LineWidth',2.5);
hold on;
plot(seg./fs, real(x(seg)),'-k','LineWidth',1.2);
hold off;
ax = gca; ax.FontSize = 24; xlabel('time(sec)')

