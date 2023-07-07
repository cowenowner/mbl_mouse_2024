lfp_bin_file_path = 'E:\NSB_Mouse\23241\01\realtest101_23241_g0\realtest101_23241_g0_imec0\realtest101_23241_g0_t0.imec0.lf.bin';
all_channels = 1:385;
[path_name, tmp] = fileparts(lfp_bin_file_path);
lfp_fname = [tmp '.bin'];
meta_fname = strrep(lfp_fname,'.lf.bin','.lf.meta');
% Extract_varargin;
obj = SGLX_Class;
LFP.meta = obj.ReadMeta(meta_fname,path_name);
LFP.nChan = str2double(LFP.meta.nSavedChans);
LFP.sFreq = str2double(LFP.meta.imSampRate);
LFP.nFileSamp = str2double(LFP.meta.fileSizeBytes) / (2 * LFP.nChan);
LFP.duration_of_recording_sec = LFP.nFileSamp/LFP.sFreq;
%%
% read the data. The last rec is typically the sync pulse.
[LFP.data_uV,LFP.meta] = obj.ReadBinVolts(0,1e6,LFP.meta ,lfp_fname,path_name); %(samp0, nSamp, meta, binName, path)
[LFP.data] = obj.ReadBin(0,1e6,LFP.meta ,lfp_fname,path_name); %(samp0, nSamp, meta, binName, path)
LFP.data_uV = LFP.data_uV*1e6;
LFP.data_uV = LFP.data_uV' -  movmean(LFP.data_uV',round(LFP.sFreq)*1.5);
LFP.data_uV = LFP.data_uV';

ABV = mean(abs(LFP.data_uV),2,'omitnan');
BAD_CH_IX =  ABV > 300 | ABV < 50;
LFP.data_uV(BAD_CH_IX,:) = nan;

LFP.t_uSec = 1e6*(linspace(start_rec, end_rec, Cols(LFP.data_uV))./LFP.sFreq)';

% Median re-reference
LFP.data_uV = LFP.data_uV - median(LFP.data_uV,1,'omitnan');
% Restrict to just the channels of interest
LFP.data_uV = LFP.data_uV(hpc_channels_ix,:);
LFP.sFreq
