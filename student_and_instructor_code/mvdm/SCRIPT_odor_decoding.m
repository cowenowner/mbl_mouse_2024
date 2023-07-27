%% depth profile for Neuropixels 2.0 data
%addpath(genpath('C:\Users\mvdmlab\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\mvdm\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));

%cd('C:\data-temp\HC_2_Neuro2_g0\preprocessed');
cd('C:\Users\mvdm\Dropbox (Dartmouth College)\projects\odor-decoding\preprocessed');
load AllSpikes;
load odor_events;

% make spikes ts -- should convert into function
nCells = length(SP);

% sort by depth along probe
depths = [SP(:).neuropixels_depth_uM];
[~, sort_idx] = sort(depths, 'descend');
SP = SP(sort_idx);

myS = ts;
for iC = nCells:-1:1

    SP(iC).t = SP(iC).t_uS * 10^-6;
    S.t{iC} = SP(iC).t;
    S.label{iC} = SP(iC).cluster_id;

end
clear SP;


%% MUA-peth
cfg_Q = []; cfg_Q.smooth = 'gauss'; cfg_Q.gausswin_sd = 0.05;
Q = MakeQfromS(cfg_Q, S);

MUA = Q;
MUA = zscore_tsd(MUA);
MUA.data = nanmean(MUA.data);

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
    out = TSDpeth(cfg_peth, MUA, cue_evt(iC));
    h(iC) = plot(out, 'Color', cols(iC), 'LineWidth', 2);
    hold on;
end

set(gca, 'FontSize', 18, 'TickDir', 'out'); box off;
legend(h, {'Cue 1', 'Cue 2', 'Cue 3'}); legend boxoff;

xlabel('Time from cue onset (s)')
ylabel('z-scored MUA')

vline(0,':');

%%
param.mua_threshold = 0.5; % MUA must exceed this z-score to count
param.twin = [0 1]; % MUA threshold crossing must occur in this window relative to cue onset

cfg = [];
cfg.method = 'raw'; cfg.minlen = 0;
cfg.threshold = param.mua_threshold;
mua_thr_iv = TSDtoIV(cfg, MUA);

%% now make IVs for each trial and keep only those that overlap with the mua threshold crossings
for iC = 1:3
    this_evt = cue_evt(iC);
    this_iv = iv(param.twin(1) + this_evt.t{1}, param.twin(2) + this_evt.t{1});
    [~, keep_idx] = IntersectIV([], this_iv, mua_thr_iv);
    
    accepted_evt(iC) = this_evt;
    accepted_evt(iC).t{1} = accepted_evt(iC).t{1}(keep_idx);
end

%% then use kept trials to build encoding model
