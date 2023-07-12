%% basic SWR detection
% assumes data has been loaded, as in plot_HC_depth_profile_2.m

% get MUA
cfg_MUA = [];
cfg_MUA.tvec = lfp_tsd.tvec'; % timebase to compute MUA on
MUA = getMUA(cfg_MUA, S);

% define channels to use
cfg_SWR = [];
cfg_SWR.ripple_ch = [49 52]; % define which channels to use for ripple-band power
cfg_SWR.control_ch = [1 4]; % define which channel(s) to use as control

%% check if those are good choices
cfg_plot = [];

all_ch = cat(2, cfg_SWR.ripple_ch, cfg_SWR.control_ch);

for iCh = 1:length(all_ch)
    this_lfp = lfp_tsd;
    this_lfp.data = this_lfp.data(all_ch(iCh), :);

    cfg_f = [];
    cfg_f.f = [1 300];
    cfg_plot.lfp(iCh) = FilterLFP(cfg_f, this_lfp);
end
MultiRaster(cfg_plot, S)

%% filter and compute power
cfg_f = [];
cfg_f.f = [140 200]; % ripple band

clear ripple_ch_power
this_lfp = lfp_tsd;
this_lfp.data = this_lfp.data(cfg_SWR.ripple_ch, :);
ripple_ch_power = FilterLFP(cfg_f, this_lfp);

% get envelope and average across channels
for iCh = 1:size(ripple_ch_power.data, 1)

    ripple_ch_power.data(iCh, :) = abs(hilbert(ripple_ch_power.data(iCh, :)));

end
ripple_ch_power.data = nanmean(ripple_ch_power.data, 1);
ripple_ch_power.data = medfilt1(ripple_ch_power.data, 11); % median filter smoothing -- could tweak this

ripple_ch_power = zscore_tsd(ripple_ch_power);


%% detect events
cfg = [];
cfg.method = 'raw';
cfg.threshold = 2;
cfg.operation =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.05; % merge events closer than this
cfg.minlen = 0.05; % minimum interval length
 
SWR_evt = TSDtoIV(cfg,ripple_ch_power);
 
%% to each event, add a field with the max z-scored power (for later selection)
cfg = [];
cfg.method = 'max'; % 'min', 'mean'
cfg.label = 'maxSWRp'; % what to call this in iv, i.e. usr.label
 
SWR_evt = AddTSDtoIV(cfg,SWR_evt,ripple_ch_power);
 
%% select only those events of >5 z-scored power
cfg = [];
cfg.operation = '>';
cfg.threshold = 3.5;
 
SWR_evt = SelectIV(cfg,SWR_evt,'maxSWRp');

%%
cfg_plot.evt = SWR_evt;
MultiRaster(cfg_plot, S)