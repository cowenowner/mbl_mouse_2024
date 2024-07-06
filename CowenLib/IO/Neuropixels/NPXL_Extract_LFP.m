function [LFP] = NPXL_Extract_LFP(lfp_bin_file_path, channels, start_end_rec )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channels are in 1 based numbering so channel 0 on Imec = channel 1 here.
% start_end_rec are also 1 based (so first record is 1, not 0).
%
% Cowen 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% start_rec = 0;
% n_samples = 100000;
% lfp_bin_file_path = fullfile(pwd,'vHP_DV7pt2_mPFC_DV7_tetrodeconfig_g0_t0.imec1.lf.bin');
if nargin <3
    start_end_rec = [];
end
if nargin <2
    channels = 1:385;
end
channels0 = channels -1;
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

TBL = NPXL_Depth_From_Meta(LFP.meta);

% grabs the middle.
end_rec = round(LFP.duration_of_recording_sec*LFP.sFreq)-1;
start_rec = round((LFP.duration_of_recording_sec-20*60)*LFP.sFreq);
samp0 = start_rec -1;

nSamp = end_rec - start_rec + 1;
LFP.start_rec = start_rec;

% read the data. The last rec is typically the sync pulse.

[LFP.data_uV] = obj.ReadBin(start_rec,nSamp,LFP.meta ,lfp_fname,path_name); %(samp0, nSamp, meta, binName, path);

nChan = str2double(LFP.meta.nSavedChans);

nFileSamp = str2double(LFP.meta.fileSizeBytes) / (2 * nChan);
samp0 = max(samp0, 0);
nSamp = min(nSamp, nFileSamp - samp0);

sizeA = [nChan, nSamp];

fid = fopen(lfp_bin_file_path, 'rb');
fseek(fid, samp0 * 2 * nChan, 'bof');
if length(channels) == 384
    % Load all
    LFP.dataArray = fread(fid, sizeA, 'int16=>double');
else
    % load a subset
    LFP.dataArray2 = nan(length(channels),nSamp,'single');

    for iCh = 1:length(channels)
        %  Go to start but also skip to the target channel.
        fseek(fid, samp0 * 2 * nChan + channels0(iCh)*2, 'bof');
        % ftell(fid)
        LFP.dataArray2(iCh,:) = fread(fid, [1 nSamp], 'int16=>double',nChan-1);
        % dataArray2(iCh,:);
        fprintf('%d,',iCh)
    end
end
fclose(fid);
clf
subplot(1,2,1)
imagesc(LFP.data_uV(:,1:20:end));colorbar
subplot(1,2,2)
imagesc(LFP.dataArray2(:,1:20:end));colorbar


% 
% LFP.data_uV   = NPXL_Convert_to_uV(LFP.data_uV,LFP.meta);
% 
% LFP.data_uV = LFP.data_uV' -  movmean(LFP.data_uV',round(LFP.sFreq)*2.5);
% LFP.data_uV = LFP.data_uV';
% 
% 
% figure; gscatter(TBL.x_inter_shank_uM,TBL.Depth_uM,TBL.shank)
% ylabel('uM'); xlabel('uM'); axis ij
% title('Shank and channel configuraiton')
% 
% 
% 
% 
% 
% 
% [path_name, tmp] = fileparts(lfp_bin_file_path);
% lfp_fname = [tmp '.bin'];
% meta_fname = strrep(lfp_fname,'.lf.bin','.lf.meta');
% % Extract_varargin;
% obj = SGLX_Class;
% LFP.meta = obj.ReadMeta(meta_fname,path_name);
% LFP.nChan = str2double(LFP.meta.nSavedChans);
% LFP.sFreq = str2double(LFP.meta.imSampRate);
% LFP.nFileSamp = str2double(LFP.meta.fileSizeBytes) / (2 * LFP.nChan);
% LFP.time_sec_of_recording = LFP.nFileSamp/LFP.sFreq;
% %LFP.meta = meta;
% LFP.start_rec = start_rec;
% % read the data. The last rec is typically the sync pulse.
% LFP.data = obj.ReadBin(start_rec,n_samples,LFP.meta ,lfp_fname,path_name); %(samp0, nSamp, meta, binName, path)
% LFP.data = int16(LFP.data ); % only a range of +/- 512 bits whcih make sense since 10bit a/d 2^10 = 1024
% % What is the conversion factor?
% v_range = 2*str2double(LFP.meta.imAiRangeMax);
% n_bits = 1024;
% mV_per_bit = v_range*1000/n_bits; % I think this ignores the gain - where is this stored?
% bits_to_mV = n_bits/(v_range*1000) % this can't be right - this has to be in the uV range, not mV
% 
% %
% figure
% subplot(2,1,1)
% imagesc(LFP.data); colorbar
% subplot(2,1,2)
% imagesc(diff(LFP.data,2,1)); colorbar % CSD  - assumes linear ordering.
% caxis([-100 100])
% % here is the reading code - will probably have to modify manually... 
% %             nChan = str2double(meta.nSavedChans);
% % 
% %             nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
% %             samp0 = max(samp0, 0);
% %             nSamp = min(nSamp, nFileSamp - samp0);
% % 
% %             sizeA = [nChan, nSamp];
% % 
% %             fid = fopen(fullfile(path, binName), 'rb');
% %             fseek(fid, samp0 * 2 * nChan, 'bof');
% %             dataArray = fread(fid, sizeA, 'int16=>double');
% %             fclose(fid);