lfp_bin_file_path = 'C:\LFPAnalysis_20230706\HC101_23242_tcat\HC101_23242_g0_tcat.imec0.lf.bin';
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

mua_filt = designfilt('bandpassiir','FilterOrder',14, 'HalfPowerFrequency1',500,'HalfPowerFrequency2',1200,'SampleRate',LFP.sFreq);
vstr_channels_ix = 0:150;
% freqz(mua_filt)
n_samples = LFP.nFileSamp;
start_rec = 0;

% read the data. The last rec is typically the sync pulse.
[LFP.data_uV,LFP.meta] = obj.ReadBinVolts(start_rec,n_samples,LFP.meta ,lfp_fname,path_name); %(samp0, nSamp, meta, binName, path)
LFP.data_uV = LFP.data_uV*1e6;
LFP.data_uV = LFP.data_uV' -  movmean(LFP.data_uV',round(LFP.sFreq)*1.5);
LFP.data_uV = LFP.data_uV';

ABV = mean(abs(LFP.data_uV),2,'omitnan');
BAD_CH_IX =  ABV > 300 | ABV < 50;
LFP.data_uV(BAD_CH_IX,:) = nan;

LFP.t_uSec = 1e6*(linspace(start_rec, LFP.nFileSamp, Cols(LFP.data_uV))./LFP.sFreq)';

% Median re-reference
LFP.data_uV = LFP.data_uV - median(LFP.data_uV,1,'omitnan');
% Restrict to just the channels of interest
LFP_vSTR = LFP.data_uV;
LFP.sFreq
