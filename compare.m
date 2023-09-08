% % % RMSE statistics of each RD methods

clear;
close all;

% load RMSE_fullyAdaptiveRD.mat
% load RMSE_RRPRD.mat
% load RMSE_MultiCurveExt.mat
load ./newARMAnoise/RMSE_MultiCurveExt.mat
load ./newARMAnoise/RMSE_RRPRD.mat
load ./newARMAnoise/RMSE_fullyAdaptiveRD.mat
load ./newARMAnoise/RMSE_SingleRD.mat

% RMSE_fullyAdaptiveRD = cell(3,3);
% for i = 1:3
%     for j = 1:3
%         RMSE_fullyAdaptiveRD{i,j} = zeros(100,1);
%     end
% end

AFUND1 = [0.1 0.2 0.5];
am = 1:3;

SNR = [Inf 5 0];

lev = 0.05 / 3 / 3 / 3;

% figure;
for i = 1:3
for snrdb = 1:3
    subplot(3,3,(i-1)*3+snrdb);
    boxchart(log10([RMSE_SingleRD{am(i),snrdb}' RMSE_fullyAdaptiveRD{am(i),snrdb} RMSE_RRPRD{am(i),snrdb} RMSE_MultiCurveExt{am(i),snrdb}]));
    xt = get(gca, 'XTick');

    pval = signrank(RMSE_SingleRD{am(i),snrdb}, RMSE_MultiCurveExt{am(i),snrdb});
    if pval<lev && median(RMSE_MultiCurveExt{am(i),snrdb})<median(RMSE_SingleRD{am(i),snrdb})
        hold on
        plot(xt(1),.6,'*k','MarkerSize',11,'LineWidth',2)
        hold off
    end

    pval = signrank(RMSE_fullyAdaptiveRD{am(i),snrdb}, RMSE_MultiCurveExt{am(i),snrdb});
    if pval<lev && median(RMSE_MultiCurveExt{am(i),snrdb})<median(RMSE_fullyAdaptiveRD{am(i),snrdb})
        hold on
        plot(xt(2),.6,'*k','MarkerSize',11,'LineWidth',2)
        hold off
    end

    pval = signrank(RMSE_RRPRD{am(i),snrdb}, RMSE_MultiCurveExt{am(i),snrdb});
    if pval<lev && median(RMSE_MultiCurveExt{am(i),snrdb})<median(RMSE_RRPRD{am(i),snrdb})
        hold on
        plot(xt(3),.6,'*k','MarkerSize',11,'LineWidth',2)
        hold off
    end
    grid on;
    ylim([-3 .6]);
%     if snrdb == 1
%         ylabel(['D = ', num2str(AFUND1(am(i)))]);
%     end
%     if i == 1
%         title(['SNR = ', num2str(SNR(snrdb))]);
%     end
    if i == 3
        set(gca, 'XTickLabel', {'\delta^{S}', '\delta^{FM}', '\delta^{RRP}', '\delta^{MH}'});
    else
        set(gca, 'XTickLabel', [])
    end
    set(gca, 'FontSize', 18);
end
end