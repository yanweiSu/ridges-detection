cd 'C:\Users\scien\Documents\TIDIS project\SST'
addpath('CGMHLK_A5data');
fff = [[1 1]]; sub = 1;
load(['CGMHLK_A5data\A5data_group_',int2str(fff(sub,1)),'.mat']);
focus = Channels{fff(sub,2),4};
label = label{fff(sub,2),1};
STAGES = cell(5,1);
for s = 1:5
    STAGES{s} = focus(:,find(label == s));
end
% input = PPG signal
PPG = STAGES{4}(:,50);
HOP = 1 ;
samplerate = 200 ;


PPGfund = [] ;
PPGmulti1 = [] ;
PPGmulti2 = [] ;
PPGmulti3 = [] ;





% start to get the true IHR (this can also be obtained by the R peaks
% of the ABP signal
LowFq = 0.01/samplerate; HighFq = 8.01/samplerate;
[~, ~, tfrsq1, ~, tfrsqtic] = ConceFT_sqSTFT_C(PPG,...
    LowFq, HighFq, 2e-4, HOP, samplerate*10+1, 1, 6, 1, 0, 0) ;
idx = find(tfrsqtic*samplerate < 3) ;


[c1] = CurveExt_M(abs(tfrsq1(idx,:)'), 0.2) ;
IHR = tfrsqtic(c1)*samplerate ;
fs2 = samplerate/HOP;

% ideally, you should see that the extracted IHR fits the dominant
% curve in the time-frequency domain.
subplot(211)
imageSQ([1:HOP:length(PPG)]/fs2, tfrsqtic*samplerate, abs(tfrsq1), .98) ;


% Do a careful observation here -- the direct multiplication of the IHR doesn't
% always fit the true IF of the multiples.
subplot(212) ; hold off
imageSQ([1:HOP:length(PPG)]/samplerate, tfrsqtic*samplerate, abs(tfrsq1), .98) ; hold on ;
plot([1:HOP:length(PPG)]/samplerate, tfrsqtic(c1)*samplerate, 'r') ;
plot([1:HOP:length(PPG)]/samplerate, tfrsqtic(c1*2)*samplerate, 'r') ;
plot([1:HOP:length(PPG)]/samplerate, tfrsqtic(c1*3)*samplerate, 'r') ;
colormap(1-gray) ;


% build fund & multiples

PPG = PPG - mean(PPG) ;

[~, ~, tfrsq1, ~, tfrsqtic] = ConceFT_sqSTFT_C(PPG,...
    LowFq, HighFq, 2e-4, HOP, samplerate*10+1, 1, 6, 1, 1, 0) ;
[h, ~] = hermf(samplerate*10+1, 1, 6) ;




c0 = round(smooth(IHR, 20, 'loess')./ ((tfrsqtic(2)-tfrsqtic(1))*samplerate)) ;

Band = round(0.4 ./ ((tfrsqtic(2)-tfrsqtic(1))*samplerate)) ;


% 1st harmonic curve extraction with the prior IHR
tmpIMG = zeros(Band*2+1, size(tfrsq1, 2)) ;
for ll = 1: size(tfrsq1, 2)
    if c0(ll)-Band >= 1
%         tmps = tfrsq1(1: c0(ll)+Band, ll) ;
%         tmps = [tmps zeros(2*Band+1-length(tmps), 1)];
%     else
        tmpIMG(:, ll) = tfrsq1(c0(ll)-Band: c0(ll)+Band, ll) ;
    end
end

[c] = CurveExt_M(abs(tmpIMG)', 0.1);
c1 = max(1, c0(:) - Band - 1) + c(:) ;
[recon] = Recon_sqSTFT_v2(tfrsq1, tfrsqtic, samplerate, c1, 0.25, h(samplerate*5+1)) ;
PPGfund = recon ;



Band = round(0.8 ./ ((tfrsqtic(2)-tfrsqtic(1))*samplerate)) ;




% 2nd harmonic curve extraction with the prior IHR
tmpIMG = zeros(Band*2+1, size(tfrsq1, 2)) ;
for ll = 1: size(tfrsq1, 2)
    if (2*c1(ll)-Band > 0 && 2*c1(ll)+Band <= size(tfrsq1, 1))
        tmpIMG(:, ll) = tfrsq1(2*c1(ll)-Band: 2*c1(ll)+Band, ll) ;
    else
        tmp = tfrsq1(max(1, 2*c1(ll)-Band): min(length(tfrsqtic), 2*c1(ll)+Band), ll) ;
        tmp = [zeros(Band*2+1 - length(tmp), 1); tmp];
    end
end

[c] = CurveExt_M(abs(tmpIMG)', 0.1);
c2 = c1(:)*2 - Band - 1 + c(:) ;
[recon] = Recon_sqSTFT_v2(tfrsq1, tfrsqtic, samplerate/HOP, c2, 0.25, h(samplerate*5+1)) ;
PPGmulti1 = recon ;




% 3rd harmonic curve extraction with the prior IHR
tmpIMG = zeros(Band*2+1, size(tfrsq1, 2)) ;
for ll = 1: size(tfrsq1, 2)
    if 3*c1(ll)-Band > 0 && 3*c1(ll)+Band <= size(tfrsq1, 1)
        tmpIMG(:, ll) = tfrsq1(3*c1(ll)-Band: 3*c1(ll)+Band, ll) ;
    else
        tmp = tfrsq1(max(1, 3*c1(ll)-Band): min(length(tfrsqtic), 3*c1(ll)+Band), ll) ;
        tmp = [zeros(Band*2+1 - length(tmp), 1); tmp] ;
    end
end

[c] = CurveExt_M(abs(tmpIMG)', 0.1);
c3 = c1(:)*3 - Band - 1 + c(:) ;
[recon] = Recon_sqSTFT_v2(tfrsq1, tfrsqtic, samplerate/HOP, c3, 0.25, h(samplerate*5+1)) ;
PPGmulti2 = recon ;

%%
if 1
    figure ;
    subplot(211)
    imageSQ([1:HOP:length(PPG)]/samplerate(1), tfrsqtic*samplerate, abs(tfrsq1), .98) ;
    axis([-inf inf 0 12]) ; set(gca,'fontsize', 20) ; ylabel('Freq (Hz)')


    subplot(212) ; hold off
    imageSQ([1:HOP:length(PPG)]/samplerate, tfrsqtic*samplerate, abs(tfrsq1), .98) ; hold on ;
    plot([1:HOP:length(PPG)]/samplerate, tfrsqtic(c1)*samplerate, 'r') ;
    plot([1:HOP:length(PPG)]/samplerate, tfrsqtic(c2)*samplerate, 'r') ;
%     plot([1:HOP:length(PPG)]/samplerate, tfrsqtic(c3)*samplerate, 'r') ;
    colormap(1-gray) ; hold off ; axis([-inf inf 0 12])
    set(gca,'fontsize', 20) ; xlabel('Time (min)') ; ylabel('Freq (Hz)')


    export_fig([num2str(PVPgoodCASEs(qqidx)) '-' num2str(qq) '-reflective-TFR'], '-transparent', '-m2');
end






%%

save([filenamePVP(1:end-4),'_fund.mat'], 'PPGfund', 'PPGmulti1',...
    'PPGmulti2', 'IHR') ;

