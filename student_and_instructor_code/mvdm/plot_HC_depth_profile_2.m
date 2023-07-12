%% depth profile for Neuropixels 2.0 data
addpath(genpath('C:\Users\mvdmlab\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));

cd('C:\data-temp\HC_2_Neuro2_g0\preprocessed');
load AllSpikes;
load odor_events;
load LFPs;

%% keep LFPs from one column of sites only
figure;
plot(lfp_tsd.usr.xcoord, lfp_tsd.usr.ycoord, '.');
title('channelmap')

% 277 309 527 559 are the x-coords
keep_idx = find(lfp_tsd.usr.xcoord == 277);
depths = lfp_tsd.usr.ycoord(keep_idx);
[sorted_depths, sort_idx] = sort(depths, 'descend');
keep_idx = keep_idx(sort_idx); % idxs of column, sorted by depth

lfp_tsd.data = lfp_tsd.data(keep_idx, :);
lfp_tsd.usr.xcoord = lfp_tsd.usr.xcoord(keep_idx);
lfp_tsd.usr.ycoord = lfp_tsd.usr.ycoord(keep_idx);

%% set LFPs up for plotting
cfg_mr = [];

cfg_f = [];
cfg_f.f = [1 300];

chan_list = 1:3:size(lfp_tsd.data, 1);
for iCh = 1:length(chan_list)

    this_lfp = lfp_tsd;
    this_lfp.data = this_lfp.data(chan_list(iCh), :);
    
    cfg_mr.lfp(iCh) = FilterLFP(cfg_f, this_lfp);

end

%% make spikes ts
nCells = length(SP);

% sort by depth along probe
depths = [SP(:).neuropixels_depth_uM];
[~, sort_idx] = sort(depths, 'descend');
SP = SP(sort_idx);

% convert to s 
myS = ts;
for iC = nCells:-1:1

    SP(iC).t = SP(iC).t_uS * 10^-6;
    S.t{iC} = SP(iC).t;
    S.label{iC} = SP(iC).cluster_id;

end

%% evt
cfg_mr.evt = [];
cfg_mr.evt = evt;
%cfg_mr.evt = SelectTS([], evt, 1);
%cfg_mr.evt{2} = SelectTS([], evt, 2);
%cfg_mr.evt{3} = SelectTS([], evt, 3);

%% plot
cfg_mr.lfpHeight = 15;
cfg_mr.lfpSpacing = 5;
cfg_mr.lfpMax = 20;
cfg_mr.lfpWidth = 0.5;
cfg_mr.spkColor = 'w';

MultiRaster(cfg_mr, S)

set(gca, 'Color', [0 0 0]);
set(gca, 'XColor', [1 1 1]);
set(gca, 'YColor', [1 1 1]);
set(gcf, 'Color', [0 0 0]);
set(gca, 'YTick', [])
ylabel('')