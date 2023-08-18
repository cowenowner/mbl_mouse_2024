%% adding cross-validation to decoding
% load the data
addpath(genpath('C:\Users\mvdmlab\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));
addpath('C:\Users\mvdmlab\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\mvdm');
%addpath(genpath('C:\Users\mvdm\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));

cd('C:\data-temp\HC_2_Neuro2_g0\preprocessed');
%cd('C:\Users\mvdm\Dropbox (Dartmouth College)\projects\odor-decoding\preprocessed');
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
    % remove illegal events
    tstart = cue_evt(iC).t{1} + cfg_peth.window(1); tend = cue_evt(iC).t{1} + cfg_peth.window(2);
    remove_idx = tstart < MUA.tvec(1) | tend > MUA.tvec(end);
    cue_evt(iC).t{1}(remove_idx) = [];

    out = TSDpeth(cfg_peth, MUA, cue_evt(iC));
    h(iC) = plot(out, 'Color', cols(iC), 'LineWidth', 2);
    hold on;
end

set(gca, 'FontSize', 18, 'TickDir', 'out'); box off;
legend(h, {'Cue 1', 'Cue 2', 'Cue 3'}); legend boxoff;

xlabel('Time from cue onset (s)')
ylabel('z-scored MUA')

vline(0,':');

%% rearrange events as table
[evt_times, evt_labels] = TSToTable(evt);
tstart = evt_times + cfg_peth.window(1); tend = evt_times + cfg_peth.window(2);
remove_idx = tstart < MUA.tvec(1) | tend > MUA.tvec(end);
evt_times(remove_idx) = []; evt_labels(remove_idx) = [];

%% construct Q-mat with event codes
cfg_Q = [];
cfg_Q.dt = 0.05;
cfg_Q.boxcar_size = 7;

peth_edges = -1:cfg_Q.dt:2; % in s, relative to cues
peth_centers = peth_edges(1:end-1)+cfg_Q.dt/2;
nBins = length(peth_edges) - 1;

big_Q = []; big_Q_bin = []; big_Q_trial = []; big_Q_cue = []; big_Q_tvec = [];

nEvt = length(evt_times);
for iEvt = 1:nEvt

    cfg_Q.tvec_edges = peth_edges + evt_times(iEvt);
    this_Q = MakeQfromS(cfg_Q, S);

    big_Q = cat(2, big_Q, this_Q.data);
    big_Q_bin = cat(2, big_Q_bin, 1:nBins);
    big_Q_trial = cat(2, big_Q_trial, iEvt.*ones(1, nBins));
    big_Q_cue = cat(2, big_Q_cue, evt_labels(iEvt).*ones(1, nBins));
    big_Q_tvec = cat(2, big_Q_tvec, cfg_Q.tvec_edges(1:end-1)+cfg_Q.dt/2);

end

big_Q_tsd = tsd;
big_Q_tsd.data = big_Q; big_Q_tsd.tvec = big_Q_tvec;
big_Q_tsd.usr.bin = big_Q_bin; big_Q_tsd.usr.trial = big_Q_trial; big_Q_tsd.usr.cue = big_Q_cue;

%% create trial test/train partitions
C = cvpartition(nEvt, 'LeaveOut');

p_map = nan(size(big_Q_tsd.tvec)); % output variable

for iP = C.NumTestSets:-1:1

    % training: construct tuning curves
    this_train_trial = find(C.training(iP));
    keep = find(ismember(big_Q_tsd.usr.trial, this_train_trial));
    this_train_Q = SelectTSD([], big_Q_tsd, keep);

    tc = QtoTC_byCue(this_train_Q);

    % testing: decode
    this_test_trial = find(C.test(iP));
    keep = find(ismember(big_Q_tsd.usr.trial, this_test_trial));
    this_test_Q = SelectTSD([], big_Q_tsd, keep);

    for iBin = 1:nBins
        tc_decode = cat(2, tc{1}(:, iBin), tc{2}(:, iBin), tc{3}(:, iBin));

        cfg_dec = [];
        p = DecodeZ(cfg_dec, this_test_Q, tc_decode);
        [~, p_map(keep(iBin))] = max(p.data(:, iBin)); % find bin with highest probability

    end % over bins

end % over partitions

%% analyze results
%iBin = 4;

clear dec_acc dec_accCue
for iBin = 1:nBins

    bins_to_compare = find(big_Q_tsd.usr.bin == iBin);
    cm = confusionmat(p_map(bins_to_compare), big_Q_tsd.usr.cue(bins_to_compare));
    cm = cm ./ repmat(sum(cm), [3 1]);

    dec_acc(iBin) = mean(diag(cm));
    for iC = 1:3
        dec_accCue{iC}(iBin) = cm(iC, iC);
    end

    %subplot(8, 8, iBin)
    %imagesc(cm); colorbar; set(gca, 'TickDir', 'out');
    %title(sprintf('%.1fs %.1f%%', peth_centers(iBin), dec_acc(iBin)*100));

end

cols = 'rgb';
figure;
for iC = 1:3
    plot(peth_centers, dec_accCue{iC}, 'Color', cols(iC), 'LineWidth', 1, 'LineStyle', '--'); hold on;
end
plot(peth_centers, dec_acc, 'k', 'LineWidth', 2); box off;
hold on;
plot(peth_centers, dec_acc, '.k', 'MarkerSize', 20);

ylim([0 1])
set(gca, 'FontSize', 18, 'TickDir', 'out')
vline(0, 'k'); 
hline(0.33, ':r');
ylabel('decoding accuracy'); xlabel('time from cue (s)')

%% now decode all data
[max_val, max_ind] = max(dec_acc) % find which time bin to use for decoding

tc = QtoTC_byCue(big_Q_tsd);
tc_decode = cat(2, tc{1}(:, max_ind), tc{2}(:, max_ind), tc{3}(:, max_ind));

cfg_Q = [];
cfg_Q.dt = 0.025;
Q_decode = MakeQfromS(cfg_Q, S);

cfg_dec = [];
p = DecodeZ(cfg_dec, Q_decode, tc_decode);

%% set up LFPs for plotting
load LFPs;

keep_idx = find(lfp_tsd.usr.xcoord == 277);
depths = lfp_tsd.usr.ycoord(keep_idx);
[sorted_depths, sort_idx] = sort(depths, 'descend');
keep_idx = keep_idx(sort_idx); % idxs of column, sorted by depth

lfp_tsd.data = lfp_tsd.data(keep_idx, :);
lfp_tsd.usr.xcoord = lfp_tsd.usr.xcoord(keep_idx);
lfp_tsd.usr.ycoord = lfp_tsd.usr.ycoord(keep_idx);

cfg_mr = [];

cfg_f = [];
cfg_f.f = [1 300];

chan_list = 1:3:size(lfp_tsd.data, 1);
for iCh = 1:length(chan_list)

    this_lfp = lfp_tsd;
    this_lfp.data = this_lfp.data(chan_list(iCh), :);
    
    cfg_mr.lfp(iCh) = FilterLFP(cfg_f, this_lfp);

end

MultiRaster(cfg_mr, S);

p1 = tsd(p.tvec, p.data(1,:));
p2 = tsd(p.tvec, p.data(2,:));
p3 = tsd(p.tvec, p.data(3,:));

ph(1) = plot(p1.tvec, p1.data*10, 'r');
ph(2) = plot(p2.tvec, p2.data*10, 'g');
ph(3) = plot(p3.tvec, p3.data*10, 'b');

%%
%%% OUTTAKES %%%

%%
param.mua_threshold = 0.5; % MUA must exceed this z-score to count
param.twin = [0 1]; % MUA threshold crossing must occur in this window relative to cue onset

cfg = [];
cfg.method = 'raw'; cfg.minlen = 0;
cfg.threshold = param.mua_threshold;
mua_thr_iv = TSDtoIV(cfg, MUA);

%% now make IVs for each trial and keep only those that overlap with the mua threshold crossings
% currently broken because relies on some version of IntersectIV() that I
% don't have here
for iC = 1:3
    this_evt = cue_evt(iC);
    this_iv = iv(param.twin(1) + this_evt.t{1}, param.twin(2) + this_evt.t{1});
    [~, keep_idx] = IntersectIV([], this_iv, mua_thr_iv);
    
    accepted_evt(iC) = this_evt;
    accepted_evt(iC).t{1} = accepted_evt(iC).t{1}(keep_idx);
end

%% then use kept trials to build encoding model
cfg_Q = []; 
cfg_Q.smooth = 'gauss'; cfg_Q.gausswin_sd = 0.1;
cfg_Q.dt = 0.1;
Q = MakeQfromS(cfg_Q, S);

nDecBins = 5; % number of decoding bins past center to use
tc = [];

cfg_peth.mode = 'raw';
for iC = 1:3
    out(iC) = TSDpeth(cfg_peth, Q, cue_evt(iC));
    outZ(iC) = zscore_tsd(out(iC));

    subplot(3,2,(iC-1)*2 + 1)
    imagesc(out(iC).tvec, 1:nCells, out(iC).data); vline(0,'k');
    set(gca, 'FontSize', 18, 'TickDir', 'out'); box off;
    xlabel(sprintf('Time from odor %d onset (s)', iC));
    title('raw')

    subplot(3,2,(iC-1)*2 + 2)
    imagesc(out(iC).tvec, 1:nCells, outZ(iC).data); vline(0,'k');
    set(gca, 'FontSize', 18, 'TickDir', 'out'); box off;
    xlabel(sprintf('Time from odor %d onset (s)', iC));
    title('z-scored')

    tc_start_bin = ceil(length(out(iC).tvec) / 2);
    %tc = cat(2, tc, out(iC).data(:, tc_start_bin:tc_start_bin + nDecBins - 1)); % assemble decoding matrix
    tc = cat(2, tc, nanmean(out(iC).data(:, tc_start_bin:tc_start_bin + nDecBins - 1), 2)); % assemble decoding matrix

end

%% decode
cfg_dec = [];
p = DecodeZ(cfg_dec, Q, tc);
[~, p_map] = max(p.data, [], 1); % find bin with highest probability

%% construct ground truth tsd
true_out = nan(size(Q.tvec));
for iC = 1:3
    for iEvt = 1:length(cue_evt(iC).t{1})

        this_idx = nearest_idx3(cue_evt(iC).t{1}(iEvt),Q.tvec);
        %true_out(this_idx:this_idx+nDecBins - 1) = (1:5) + (iC-1)*nDecBins;
        true_out(this_idx) = iC;

    end
end

cm = confusionmat(p_map, true_out);