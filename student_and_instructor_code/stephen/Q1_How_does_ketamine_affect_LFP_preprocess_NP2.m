% Extracts LFPs and aligns data.
% These are neuropixels 2 data.
%  C:\Data\NSandB_course2023\23239\01_7_10_23\vSTR_Ket_Neuro2_g0
%% DESPITE THE NAME, this is NOT neuropixels 2.0 - it's 1.0
clearvars
decimation_factor = 60; % 30000/120 = 250;
channel_skip = 6;

root_data_dir = 'C:\Data\NSandB_course2023\23240\1_7_10_23\vSTR_Ket_Neuro2_g0';
[~,n] = fileparts(root_data_dir);
ap_data_dir = fullfile(root_data_dir,[n '_imec0']);
% Get the channel map so that we know the depths of each site.
cd (ap_data_dir)
% SGLXMetaToCoords();
% Load the data.
d = dir(fullfile(ap_data_dir,'*.ap.bin'));
lfp_bin_file_path = fullfile(ap_data_dir,d(1).name);
[path_name, tmp] = fileparts(lfp_bin_file_path);
lfp_fname = [tmp '.bin'];
meta_fname = strrep(lfp_fname,'.ap.bin','.ap.meta');
obj = SGLX_Class;
LFP.meta = obj.ReadMeta(meta_fname,path_name);
LFP.nChan = str2double(LFP.meta.nSavedChans);
LFP.sFreq = str2double(LFP.meta.imSampRate);
LFP.nFileSamp = str2double(LFP.meta.fileSizeBytes) / (2 * LFP.nChan);
LFP.time_sec_of_recording = LFP.nFileSamp/LFP.sFreq;
start_rec = 0;
n_samples = LFP.nFileSamp;
% read the data. The last rec is typically the sync pulse.
% read the data in blocks and downsample.
% Fast downsample and save the bin file for processing.
NPXL_Downsample_And_Save_LFP(lfp_fname, LFP.meta , decimation_factor, channel_skip);