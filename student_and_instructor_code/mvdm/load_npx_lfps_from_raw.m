%%
addpath('C:\Users\mvdmlab\Documents\GitHub\neuropixel-utils');

%setenv('NPIX_MAP_FILE', 'C:\Users\mvdmlab\Documents\GitHub\neuropixel-utils\map_files\neuropixPhase3B2_kilosortChanMap.mat') % wrong but works?
setenv('NPIX_MAP_FILE', 'C:\data-temp\HC_2_Neuro2_g0\HC_2_Neuro2_g0_imec0\HC_2_Neuro2_g0_tcat.imec0.ap_kilosortChanMap.mat');
imec = Neuropixel.ImecDataset(pwd);
mmap = imec.memmapAP_full();

load('HC_2_Neuro2_g0_tcat.imec0.ap_kilosortChanMap.mat');

%% load LFPs into MultiRaster format
Fs = 30000;

cfg_f = [];
cfg_f.f = [1 200];

ch_list = 20:20:360;
for iCh = 1:length(ch_list)

    fprintf('%d/%d...\n',iCh,length(ch_list));

    this_lfp = mmap.Data.x(ch_list(iCh), :);
    tvec = 1/Fs * (0:1:length(this_lfp)-1);

    this_lfp = decimate(double(this_lfp), 30);
    tvec = tvec(1:30:end);

    this_tsd = tsd(tvec,this_lfp);

    this_tsd.cfg.hdr{1}.Fs = Fs / 30;

    cfg_mr.lfp(iCh) = FilterLFP(cfg_f,this_tsd);

end

%% load LFPs into tsd for general use
Fs = 30000;

cfg_f = [];
cfg_f.f = [1 475];
decimate_factor = 15;

nSamples = size(mmap.Data.x,2);
tvec = 1/Fs * (0:1:nSamples-1);
tvec = tvec(1:decimate_factor:end);

ch_list = 1:384;

all_lfps = zeros(length(ch_list), length(tvec));

for iCh = length(ch_list):-1:1

    fprintf('%d/%d...\n',iCh,length(ch_list));

    this_lfp = mmap.Data.x(ch_list(iCh), :);
    all_lfps(iCh, :) = decimate(double(this_lfp), decimate_factor);
    
    usr.xcoord(ch_list(iCh)) = xcoords(ch_list(iCh));
    usr.ycoord(ch_list(iCh)) = ycoords(ch_list(iCh));
    
end

this_tsd = tsd(tvec,all_lfps);
this_tsd.cfg.hdr{1}.Fs = Fs / decimate_factor;
this_tsd.usr = usr;
    