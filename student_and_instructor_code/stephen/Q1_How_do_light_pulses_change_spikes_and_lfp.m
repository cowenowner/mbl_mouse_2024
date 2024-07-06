% Do LED blinky lights alter brain activity in the hippocampus?
%% 
clearvars
data_folder = 'D:\Data\';
npxl_top_dir_name = 'PhotoPixelsStrobe_g0';
kilosort_dir_name = 'kilosort_cowen';

% Figure out directories and filenames
[D] = WH_load_NPXL_data(data_folder,npxl_top_dir_name,kilosort_dir_name );

% Load data
LFP = NPXL_Extract_LFP (D.lfp_bin_file_path,1:20:385);
obj = SGLX_Class;

LFP.meta = obj.ReadMeta(D.meta_fname_lf, D.ap_data_dir );
LFP.nChan = str2double(LFP.meta.nSavedChans);
LFP.sFreq = str2double(LFP.meta.imSampRate);
LFP.nFileSamp = str2double(LFP.meta.fileSizeBytes) / (2 * LFP.nChan);
LFP.time_sec_of_recording = LFP.nFileSamp/LFP.sFreq;
% read the data. The last rec is typically the sync pulse.
% read the data in blocks and downsample.
% Fast downsample and save the bin file for processing.
NPXL_Downsample_And_Save_LFP(D.lfp_bin_file_path, LFP.meta , 2, 4);
