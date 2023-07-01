function TI = DANA_trial_info_from_stim_times(stim_times_sec)
% Extract the trial start and end and LV for each trial condition.
% Cowen 2022
inter_trial_time_sec = 100; % reasonable assumption
d = diff(stim_times_sec);
ix = find(d>inter_trial_time_sec);
TI.end_times_sec = stim_times_sec(ix);
TI.end_times_sec(end+1) = stim_times_sec(end);
tmp_sec = [0; stim_times_sec(:)];
d2 = diff(tmp_sec);
ix2 = find(d2>inter_trial_time_sec);
TI.start_times_sec = stim_times_sec(ix2);
% For each trial, calculate the local variance
for iT = 1:length(TI.start_times_sec)
    TI.TrialID(iT) = iT;
    IX = stim_times_sec >= TI.start_times_sec(iT) & stim_times_sec <= TI.end_times_sec(iT);
    t = stim_times_sec(IX);
    [TI.LV(iT),TI.LVR(iT)] = LocalVariance(diff(t)*1000);
    TI.Hz(iT) = length(t)/(t(end)-t(1));
    TI.Stim_sequence_sec{iT} = t - TI.start_times_sec(iT);
end
%
TI.TrialID = TI.TrialID(:);
TI.LV = TI.LV(:);
TI.LVR = TI.LVR(:);
TI.Hz = TI.Hz(:);

% figure out a within-condition trial ID 
TI.within_condition_trial_ID = zeros(size(TI.Hz));
TI.Hz_group = nan(size(TI.Hz));

for iHz = [5 10 20 40 60]
    TI.Hz_group(TI.Hz > iHz-2 & TI.Hz < iHz + 2) = iHz;
end

u = unique(TI.Hz_group);
for ii = 1:length(u)
    ix = find(TI.Hz_group == u(ii));
    TI.within_condition_trial_ID(ix) = 1:length(ix);
end

TI.LV_group = nan(size(TI.Hz));
for iLV = [0 .3 1 1.2]
    TI.LV_group(TI.LV > iLV -.1 & TI.LV < iLV + .1) = iLV;
end


if nargout == 0
    figure
    plot(stim_times_sec,ones(size(stim_times_sec)),'.')
    hold on
    plot(TI.start_times_sec,ones(size(TI.start_times_sec)),'g>')
    plot(TI.end_times_sec,ones(size(TI.end_times_sec)),'r<')

    figure
    plot_raster(TI.Stim_sequence_sec)
    xlabel('sec')
    ylabel('trial')
    axis ij
    pubify_figure_axis
    title('Stim sequences')
end