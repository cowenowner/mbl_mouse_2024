%% load photometry data
cd('E:\Ach_DualPixel_7_14_23_g0\Ach_DualPixel_7_14_23_g0_imec1')
restoredefaultpath;
addpath(genpath('C:\Users\mvdmlab\Documents\GitHub\mbl_mouse_2023\CowenLib'));

obj = SGLX_Class;

[binName,path] = uigetfile('*.bin', 'Select Binary File');

meta = SGLX_Class.ReadMeta(binName, path);
dataArray = SGLX_Class.ReadBin(0, Inf, meta, binName, path);

this_data = dataArray(8, :);
dec_factor = 25;
this_data = decimate(this_data, dec_factor);

tvec = 0:length(this_data)-1;
Fs = (str2num(meta.niSampRate) ./ dec_factor);
tvec = tvec .* (1 ./ Fs);

restoredefaultpath;
addpath(genpath('C:\Users\mvdmlab\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));

this_tsd = tsd(tvec,this_data);
this_tsd.cfg.hdr{1}.Fs = Fs;

cfg_f = [];
cfg_f.f = [1 20];
this_tsdf = FilterLFP(cfg_f,this_tsd);

%% add 2.0 probe data
load AllSpikes;
nCells = length(SP);

load('Ach_DualPixel_7_14_23_g0_tcat.imec1.ap_kilosortChanMap.mat')

% sort by depth along probe
% remember that MultiRaster plots cell 1 at the bottom, so if depth is sorted in descending order, then deepest cell is cell 1
depths = [SP(:).neuropixels_depth_uM];
[~, sort_idx] = sort(depths, 'descend'); 
SP = SP(sort_idx);

% now separate by shank
sort_idx = [];
%cfg_mr = [];
cfg_mr.spkColor = zeros(nCells, 3);
cols = [0 0 0; 0.7 0 0; 0 0.7 0; 0 0 0.7];
for iShank = 1:4

    this_shank_chan = chanMap0ind(kcoords == iShank);

    % now need to know which SP.template_index0 are in this set
    this_idx = find(ismember([SP.template_index0], this_shank_chan));
    sort_idx = cat(2, sort_idx, this_idx);
    cfg_mr.spkColor(this_idx, :) = repmat(cols(iShank, :), [length(this_idx) 1]);

end
SP = SP(sort_idx);
cfg_mr.spkColor = cfg_mr.spkColor(sort_idx, :);

% convert to s 
S = ts;
for iC = nCells:-1:1

    SP(iC).t = SP(iC).t_uS * 10^-6;
    S.t{iC} = SP(iC).t;
    S.label{iC} = SP(iC).cluster_id;

end

cfg_mua = [];
cfg_mua.tvec = this_tsdf.tvec;
cfg_mua.sigma = 0.01;
this_mua = getMUA(cfg_mua,S);
this_mua.data = this_mua.data';

%% load LFPs

% find top and bottom channel of each shank
ymin = min(ycoords);
ymax = max(ycoords);
lfp_idx = []; lfp_col = [];
for iShank = 4:-1:1 % for LFPs, plotting order is opposite compared to spikes: each LFP is plotted below the previous one
    lfp_idx = cat(2, lfp_idx, find(kcoords == iShank & ycoords == ymax));
    lfp_idx = cat(2, lfp_idx, find(kcoords == iShank & ycoords == ymin));
    lfp_col = cat(2, lfp_col, [iShank iShank]);
end
%lfp_idx = chanMap0ind(lfp_idx(1,:)); % two sites per shank are equal top and bottom, so pick one
lfp_idx = lfp_idx(1,:); % hmm when loading it seems there isn't a channel 0 -- that's scary

addpath('C:\Users\mvdmlab\Documents\GitHub\neuropixel-utils');

setenv('NPIX_MAP_FILE', 'E:\Ach_DualPixel_7_14_23_g0\Ach_DualPixel_7_14_23_g0_imec1\Ach_DualPixel_7_14_23_g0_tcat.imec1.ap_kilosortChanMap.mat');
imec = Neuropixel.ImecDataset(pwd);
mmap = imec.memmapAP_full();

Fs = 30000;

cfg_f = [];
cfg_f.f = [1 200];

for iCh = 1:length(lfp_idx)

    fprintf('%d/%d...\n',iCh,length(lfp_idx));

    this_lfp = mmap.Data.x(lfp_idx(iCh), :);
    tvec = 1/Fs * (0:1:length(this_lfp)-1);

    this_lfp = decimate(double(this_lfp), 30);
    tvec = tvec(1:30:end);

    this_tsd = tsd(tvec,this_lfp);

    this_tsd.cfg.hdr{1}.Fs = Fs / 30;

    cfg_mr.lfp(iCh + 2) = FilterLFP(cfg_f,this_tsd);

end

%% now add 1.0 probe data
cd('E:\Ach_DualPixel_7_14_23_g0\Ach_DualPixel_7_14_23_g0_imec0')

load AllSpikes;
nCells2 = length(SP);

depths = [SP(:).neuropixels_depth_uM];
[~, sort_idx] = sort(depths, 'descend'); 
SP = SP(sort_idx);

SP = SP(sort_idx);

% convert to s 
%for iC = nCells+1:1:nCells+nCells2
for iC = 1:nCells2

    SP(iC).t = SP(iC).t_uS * 10^-6;
    S.t{nCells+iC} = SP(iC).t;
    S.label{nCells+iC} = SP(iC).cluster_id;
    cfg_mr.spkColor(nCells+iC, :) = [0.7 0.7 0.7];

end
%%
%cfg_mr = [];
cfg_mr.lfp(1) = this_mua; cfg_mr.lfpColor(1, :) = [0.7 0.7 0.7];
cfg_mr.lfp(2) = this_tsdf; cfg_mr.lfpColor(2, :) = [0 1 1];
cfg_mr.lfpColor(3:3+length(lfp_idx)-1, :) = cols(lfp_col, :);

MultiRaster(cfg_mr, S);

%% cell x-corrs with Ach

% first make each cell's tsd firing rate, then xcorr with ACh

all_xc = [];
for iC = nCells:-1:1

    this_S = SelectTS([], S, iC);

    cfg_mua = [];
    cfg_mua.tvec = cfg_mr.lfp(2).tvec;
    cfg_mua.sigma = 0.05;
    this_mua = getMUA(cfg_mua, this_S); this_mua.data = this_mua.data';

    this_lfp = cfg_mr.lfp(2);

    [this_xc, xc_t] = xcorr(this_mua.data, this_lfp.data(2:end-1), 2000, 'coeff');
    all_xc(iC, :) = this_xc;

end
xc_t = xc_t ./ this_tsd.cfg.hdr{1}.Fs;

imagesc(xc_t, 1:333, all_xc)
h = vline(0)
set(gca, 'FontSize', 18, 'TickDir', 'out');
xlabel('spike timing relative to striatal ACh (s)');
ylabel('cell# in thalamus')
colorbar

%% LFP x-corrs with ACh

all_lfp_xc = [];
for iLFP = 8:-1:1

    this_lfp = cfg_mr.lfp(iLFP+2).data;
    this_lfp_interp = interp1(cfg_mr.lfp(iLFP+2).tvec, cfg_mr.lfp(iLFP+2).data, cfg_mr.lfp(2).tvec);

    this_ach = cfg_mr.lfp(2).data;
    keep = ~isnan(this_lfp_interp);
    this_lfp_interp = this_lfp_interp(keep);
    this_ach = this_ach(keep);

    [this_xc, xc_t] = xcorr(this_lfp_interp, this_ach, 2000, 'coeff');
    all_lfp_xc(iLFP, :) = this_xc;

end

xc_t = xc_t ./ this_tsd.cfg.hdr{1}.Fs;

imagesc(xc_t, 1:8, all_lfp_xc)
h = vline(0)
set(gca, 'FontSize', 18, 'TickDir', 'out');
xlabel('LFP voltage relative to striatal ACh (s)');
ylabel('LFP# in thalamus')
colorbar