%%
% script to extract odor presentation events
% requires *g0_tcat.nidq.xia_x_0.txt file with TTL times
% requires vandermeerlab codebase

%% params
cfg = [];
cfg.line_noise_isi = [1/70 1/50]; % remove any ISIs in this range
cfg.true_isi = [4.5 5.5]; % ISIs in this range are true trials
cfg.code_isi = [0.05 0.1]; % odor identity coding pulses

%% load
fn = FindFile('scowen*.txt'); % TTL pulse times that Stephen extracted with his magical CatGT command
evt_t = textread(fn);

% remove line noise
diffs = diff(evt_t);
line_noise_diff_idx = find(diffs > cfg.line_noise_isi(1) & diffs <= cfg.line_noise_isi(2));
fprintf('Removing %d/%d events due to line noise...\n', length(line_noise_diff_idx)+1, length(evt_t));

remove_idx = [line_noise_diff_idx-1; line_noise_diff_idx];
evt_t(remove_idx) = [];

% if two diffs are wrong in a row, remove data point corresponding to first
% of those diffs
diffs = diff(evt_t);
keep = (diffs > cfg.true_isi(1) & diffs <= cfg.true_isi(2)) | (diffs > cfg.code_isi(1) & diffs <= cfg.code_isi(2));
    
prev_keep = 1;
remove_idx = [];
for iK = 1:length(keep)

    this_keep = keep(iK);
    if this_keep == 0 & prev_keep == 0
        remove_idx = [remove_idx iK];
    end
    prev_keep = this_keep;
    
end 
evt_t(remove_idx) = [];   
fprintf('Removed %d rogue events...\n', length(remove_idx));

% next, remove unknown events at edges iteratively
totalRemoved = 0;
iRemoved = 1;
while iRemoved > 0
    
    diffs = diff(evt_t);
    
    % these diffs we think are important
    keep = (diffs > cfg.true_isi(1) & diffs <= cfg.true_isi(2)) | (diffs > cfg.code_isi(1) & diffs <= cfg.code_isi(2));
    
    iRemoved = 0;
    
    if keep(1) == 0 % first diff wrong
        iRemoved = 1;
        evt_t = evt_t(2:end);
    end
    
    if keep(end) == 0 % last diff wrong
        iRemoved = 1;
        evt_t = evt_t(1:end-1);
    end
    
    totalRemoved = totalRemoved + iRemoved;
    
end
fprintf('Removed %d unknown edge events...\n', totalRemoved);

%% now classify events according to number of pulses
diffs = diff([evt_t(1)-5; evt_t]);
trial_diff_idx = [find(diffs > cfg.true_isi(1) & diffs <= cfg.true_isi(2))];

trial_times = evt_t(trial_diff_idx);

pulse_numbers = diff([trial_diff_idx; length(evt_t)+1]); 

%%
evt = ts;
for iT = 1:3
   evt.t{iT} = trial_times(pulse_numbers == iT);
   evt.label{iT} = num2str(iT);
end