function [cALL] = MultiCurveExt(tfr, tfrtic, lambda, mu, num_multi, band_fund, band_multi, version, tt, bw)

tfrtic = reshape(tfrtic,[],1);
cALL = [];
flag = 0;

if version==2 && length(lambda)==2 && length(mu)==1 && length(num_multi)==1 && length(band_multi)==1

fprintf("Run 2-curves simultaneous extraction\n")
% mex './CurveMultiExt_init_2curves.c';
c_init = [1; 1];
for t = 1:length(tt)-1
    TFR = abs(tfr(:, tt(t)+1:tt(t+1)));
    c0 = ones(size(TFR,2),1);
    c1 = ones(size(c0));
    fprintf("%d/%d: %.2f Hz ~ %.2f Hz:\n", t, length(tt)-1, band_fund(1), band_fund(2));

    tic
    [c0, c1] = CurveMultiExt_init_2curves(TFR.', tfrtic, lambda, mu, num_multi, band_fund, band_multi, flag, c_init);
    toc

    c_init = [c0(end); c1(end)];
    band_fund = [max(tfrtic(c0(end))+bw(1), 0), min(tfrtic(c0(end))+bw(2), tfrtic(end))];
    cALL = [cALL; [c0 c1]];
    flag = 1;
end

elseif version == 3 && length(lambda)==3 && length(mu)==2 && length(num_multi)==2 && length(band_multi)==2

fprintf("Run 3-curves simultaneous extraction\n")
% mex './CurveMultiExt_init_3curves.c';
c_init = [1; 1; 1];
for t = 1:length(tt)-1
    TFR = abs(tfr(:, tt(t)+1:tt(t+1)));
    c0 = ones(size(TFR,2),1);
    c1 = ones(size(c0));
    c2 = ones(size(c0));
    fprintf("%d/%d: %.2f Hz ~ %.2f Hz:\n", t, length(tt)-1, band_fund(1), band_fund(2));
    
    tic
    [c0, c1, c2] = CurveMultiExt_init_3curves(TFR.', tfrtic, lambda, mu, num_multi, band_fund, band_multi, flag, c_init);
    toc

    c_init = [c0(end); c1(end); c2(end)];
    band_fund = [max(tfrtic(c0(end))+bw(1), 0), min(tfrtic(c0(end))+bw(2), tfrtic(end))];
    cALL = [cALL; [c0 c1 c2]];
    flag = 1;
end

else

fprintf("No such thing!\n");
return;

end

end