%% multi-unit activity (sum of all spikes over time)
cfg_MUA = [];
cfg_MUA.tvec = lfp_tsd.tvec'; % timebase to compute MUA on
MUA = getMUA(cfg_MUA, S);

figure; 
subplot(121)
plot(MUA)
title('raw MUA (spks/s)')

subplot(122)
MUAz = zscore_tsd(MUA);
plot(MUAz)
title('z-scored MUA (a.u.)');


%% event-triggered MUA
for iC = 1:3
    cue_evt(iC) = SelectTS([], evt, iC); % select events for each individual cue
end

cfg_peth = []; % parameters for PETH
cfg_peth.window = [-1 1];
cfg_peth.dt = 0.01;
cfg_peth.mode = 'interp';

figure;

cols = 'rgb';
for iC = 1:3
    out = TSDpeth(cfg_peth, MUAz, cue_evt(iC));
    h(iC) = plot(out, 'Color', cols(iC), 'LineWidth', 2);
    hold on;
end

set(gca, 'FontSize', 18, 'TickDir', 'out'); box off;
legend(h, {'Cue 1', 'Cue 2', 'Cue 3'}); legend boxoff;

xlabel('Time from cue onset (s)')
ylabel('z-scored MUA')

vline(0,':');

%% some single-unit PETHs
for iC = 1:9

    this_S = SelectTS([], S, iC); % select only the current cell
    MUA = getMUA(cfg_MUA, this_S); % "MUA" for one cell is just that cell's firing rate
    MUAz = zscore_tsd(MUA);

    subplot(3,3,iC)

    for iCue = 1:3
        this_out = TSDpeth(cfg_peth, MUAz, cue_evt(iCue)); hold on;
        plot(this_out, 'Color', cols(iCue), 'LineWidth', 2);
        title(sprintf('cell %d', iC));
    end
    vline(0,':');
end