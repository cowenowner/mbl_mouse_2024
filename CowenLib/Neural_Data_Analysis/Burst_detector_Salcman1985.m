function [BRST] = Burst_detector_Salcman1985(t_sec)
% INPUT:
%   A series of timestamps in seconds.
% OUTPUT:
%
% See https://pubmed.ncbi.nlm.nih.gov/3998798/
% UNDEVELOPED _ JUST SKETCHING OUT PIECES - NOT WORKING.
%
% 1/mean(diff(t_sec))
% Cowen 2023
% clearvars
BRST = [];
if nargin == 0
    ISI = [gamrnd(1,2,1500,1);  gamrnd(3,1,1000,1)]/15;
    ISI = ISI(randperm(length(ISI)));
    ISI(ISI<.002) = [];
    clf;histogram(log10(ISI))
    t_sec = cumsum(ISI);
    1/mean(diff(t_sec))
    figure(2); clf;plot_raster(t_sec)
end
min_consecutive_spikes = 3;

%
dur_s = t_sec(end) - t_sec(1);
mean_rate = length(t_sec)/dur_s;
initial_thresh = 1/mean_rate/2;

T = dur_s;
lambda1 = length(t_sec)/T;
n = length(t_sec);
P = poisscdf(n,lambda1*T,'upper');
S = -1*log10(P); % or is it log 10?
% Now go through and detect.
[burst_start_s, durations] = ...
    Burst_detector(t_sec,initial_thresh,min_consecutive_spikes);
BIX = durations < 0.006;
burst_start_s = burst_start_s(~BIX);
durations = durations(~BIX);
burst_end_s = burst_start_s + durations;

% Now refine using S.
figure(2);clf; plot_raster(t_sec); hold on; plot(burst_start_s,ones(size(burst_start_s)),'g>')
plot(burst_end_s,ones(size(burst_end_s)),'r<')

Sb = [];
for iB = 1:length(burst_start_s)
    dur_s = burst_end_s(iB) -  burst_start_s(iB);
    ix_st = find(t_sec == burst_start_s(iB));
    ix_ed = find(t_sec == burst_end_s(iB));
    n = ix_ed - ix_st + 1;
    Pb = poisscdf(n,lambda1*dur_s,'upper');
    aPb(iB) = Pb;
    Sb(iB) = -1*log10(Pb);
    % next step is to keep widening the window (going to the next and
    % previous spikes) until the surprise falls below a target value.
end

figure;plot(Sb);figure;plot(aPb)
