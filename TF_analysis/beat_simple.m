function beats = beat_simple(onset, osr, tempo, alpha)
% beats = beat_simple(onset, osr, tempo, alpha)
% Core of the DP-based beat tracker
% <onset> is the onset strength envelope at frame rate <osr>
% <tempo> is the target tempo (in BPM)
% <alpha> is weight applied to transition cost
% <beats> returns the chosen beat sample times (in sec).
% start: the lower bound of tempo
% 2007¡V06¡V19 Dan Ellis dpwe@ee.columbia.edu
% 2016-08-13 Revised by Li Su

onset = [zeros(1,2000) onset];
tempo = [ones(1,2000) tempo];
if nargin < 4; alpha = 100; end
% backlink(time) is best predecessor for this point
% cumscore(time) is total cumulated score to this point
localscore = onset;
backlink = -ones(1,length(localscore));
cumscore = zeros(1,length(localscore));

% period = (1/start)*osr;
% prange = round(-2*period):-round(period/2);
for i = 2001:length(localscore)%-100
    % convert bpm to samples
    period = (1/tempo(i))*osr;
    % Search range for previous beat
    prange = round(-2*period):-round(period/2);
    % Log-gaussian window over that range
    txwt = (-alpha*abs((log(prange/-period)).^2));

    timerange = i + prange;
%     if timerange(1)<1
%         timerange = timerange(timerange>0);
%     end

    % Search over all possible predecessors
    % and apply transition weighting
    scorecands = txwt + cumscore(timerange);
    % Find best predecessor beat
    [vv,xx] = max(scorecands);
    % Add on local score
    cumscore(i) = vv + localscore(i);
    % Store backtrace
    backlink(i) = timerange(xx);
end
% Start backtrace from best cumulated score
[vv,beats] = max(cumscore);

% .. then find all its predecessors
while backlink(beats(1)) > 0
    beats = [backlink(beats(1)),beats];
end
beats = beats - 2000;
beats = beats(beats>0);
% convert to seconds
% beats = (beats-1)/osr;