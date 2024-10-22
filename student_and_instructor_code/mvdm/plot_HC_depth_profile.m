%% load the data
% NOTE: this requires cowenlablibcode
%
% git clone https://github.com/cowenowner/mbl_mouse_2023
%
% then set a clean path (I use a MATLAB Favorite):
%
% restoredefaultpath;
% addpath(genpath('C:\Users\mvdmlab\Documents\GitHub\mbl_mouse_2023\CowenLib'));
lfp_bin_file_path = 'C:\Users\Administrator\Desktop\preprocessed\HC101_23242_g0_tcat.imec0.lf.bin';
all_channels = 1:385; nCh = all_channels(end);

[path_name, tmp] = fileparts(lfp_bin_file_path);
lfp_fname = [tmp '.bin'];
meta_fname = strrep(lfp_fname,'.lf.bin','.lf.meta');

epoch_time = 1800; % in s, from end of session

obj = SGLX_Class;
LFP.meta = obj.ReadMeta(meta_fname,path_name);
LFP.nChan = str2double(LFP.meta.nSavedChans);
LFP.sFreq = str2double(LFP.meta.imSampRate);
LFP.nFileSamp = str2double(LFP.meta.fileSizeBytes) / (2 * LFP.nChan);
LFP.duration_of_recording_sec = LFP.nFileSamp/LFP.sFreq;

% what data to get
end_rec = round(LFP.duration_of_recording_sec*LFP.sFreq)-1;
start_rec = end_rec - round(epoch_time*LFP.sFreq);

n_samples = end_rec - start_rec + 1;
LFP.start_rec = start_rec;

% read the data. The last rec is typically the sync pulse.
[LFP.data_uV, LFP.meta] = obj.ReadBinVolts(start_rec,n_samples,LFP.meta ,lfp_fname,path_name);

%% decimate LFPs
cfg.decimate_factor = 3;
for iCh = all_channels(end):-1:1
    LFP.data(iCh, :) = decimate(double(LFP.data_uV(iCh,:)), cfg.decimate_factor);
end
LFP.sFreq = LFP.sFreq / cfg.decimate_factor;

% create new vector of timestamps for each LFP sample
LFP.tvec = 0:size(LFP.data, 2)-1; LFP.tvec = LFP.tvec / LFP.sFreq;
LFP.tvec = -LFP.tvec + LFP.duration_of_recording_sec;
LFP.tvec = LFP.tvec(end:-1:1);

%% load spikes
cd(path_name);
load AllSpikes.mat % output file containing spike times, saved by spike sorting workflow
nCells = length(SP);

% convert to s and align with LFPs
for iC = 1:nCells

    SP(iC).t = SP(iC).t_uS * 10^-6;
    SP(iC).t = SP(iC).t(SP(iC).t > (LFP.duration_of_recording_sec - epoch_time));

end

% sort by depth along probe
depths = [SP(:).neuropixels_depth_uM];
[~, sort_idx] = sort(depths, 'descend');
SP = SP(sort_idx);

%% load events (??)

evt_fn = 'C:\data-temp\Neuropixels_HC_23242\HC101_23242_g0\synced_HC101_23242_g0_tcat.nidq.xa_6_0.txt';
evt_t = textread(evt_fn);
evt_t = evt_t(evt_t > (LFP.duration_of_recording_sec - epoch_time));

%% make mvdmlab data structures (change path)
% NOTE: this requires the vandermeerlab codebase
%
restoredefaultpath;
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));

% first choose channels to use -- these are good for our first hippcampal
% recording with the manually presented odors
cfg_plot.chan_list = 43:5:350;

% make time-stamped data (tsd) for LFPs
myLFP = tsd;
myLFP.tvec = LFP.tvec; %%%%%%%%ERROR HERE
myLFP.data = LFP.data(cfg_plot.chan_list,:);
for iLFP = 1:size(myLFP.data, 1)
    myLFP.cfg.hdr{iLFP}.Fs = LFP.sFreq;
end

% filter in LFP range, this makes the plots easier to read
cfg_f = [];
cfg_f.f = [1 200];
myLFPf = FilterLFP(cfg_f,myLFP);

% make timestamp struct for spikes
myS = ts;
for iC = nCells:-1:1

    S.t{iC} = SP(iC).t;
    S.label{iC} = SP(iC).cluster_id;

end

% make events struct
evt_ts = ts;
evt_ts.t{1} = evt_t;
evt_ts.label{1} = 'TTL';

%% set up plot
cfg_mr = [];
cfg_mr.lfpHeight = 20;
cfg_mr.lfpSpacing = 4;
cfg_mr.lfpWidth = 0.5;
cfg_mr.lfpColor = 'k';
cfg_mr.evt = evt_ts;

% arrange LFPs to plot in the way MultiRaster() expects
numLFP = size(myLFPf.data, 1);
for iLFP = 1:numLFP
    this_LFP = myLFPf;
    this_LFP.data = this_LFP.data(iLFP, :);
    cfg_mr.lfp(numLFP+1 - iLFP) = this_LFP; % reverse order, because MultiRaster plots each *subsequent* LFP below the previous one
end

% set up color for each single units
cmap = linspecer(numLFP); % RGB values for each LFP channel
LFP_chan_to_depth_on_probe = (1:384) * 10; % this is a guess -- should be stored somewhere?
these_LFP_depths = LFP_chan_to_depth_on_probe(cfg_plot.chan_list(end:-1:1)); % reverse order -- this should be done in a more controlled way
for iC = nCells:-1:1

    this_depth = SP(iC).depth_on_probe;
    [~, matching_LFP_id] = min(abs(this_depth - these_LFP_depths));
    cfg_mr.spkColor(iC, :) = cmap(matching_LFP_id , :);
end

% finally, plot
MultiRaster(cfg_mr, S);

set(gca, 'Color', [0 0 0]);
set(gca, 'XColor', [1 1 1]);
set(gca, 'YColor', [1 1 1]);
set(gcf, 'Color', [0 0 0]);
set(gca, 'YTick', [])
ylabel('')