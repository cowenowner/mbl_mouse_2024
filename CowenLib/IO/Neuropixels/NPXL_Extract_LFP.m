function [LFP] = NPXL_Extract_LFP(lfp_bin_file_path, channels, start_end_rec )
% function [LFP] = NPXL_Extract_LFP(lfp_bin_file_path, channels, start_end_rec )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lfp_bin_file_path = name of the .lf.
% channels = Channels are in 1 based numbering so channel 0 on Imec = channel 1 here.
% start_end_rec = start_end_rec are also 1 based (so first record is 1, not 0).
%
% Outputs LFP structure with the data converted to uV
%
% Cowen 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
PLOT_IT = false;
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
LFP.channels = channels;
LFP.nFileSamp = str2double(LFP.meta.fileSizeBytes) / (2 * LFP.nChan);
LFP.duration_of_recording_sec = LFP.nFileSamp/LFP.sFreq;

if nargin <3
    start_end_rec = [0 LFP.nFileSamp-LFP.nChan*2];
end
% grabs the middle.
if 0
    start_end_rec(2) = round(LFP.duration_of_recording_sec*LFP.sFreq)-1;
    start_end_rec(1) = round((LFP.duration_of_recording_sec-20*60)*LFP.sFreq);
end

TBL = NPXL_Depth_From_Meta(LFP.meta);

samp0 = start_end_rec(1) -1;
samp0 = max(samp0, 0);

LFP.nSamp = start_end_rec(2) - start_end_rec(1) + 1;
LFP.start_rec = start_end_rec(1);

% read the data. The last rec is typically the sync pulse.

% [LFP.data_uV] = obj.ReadBin(start_rec,nSamp,LFP.meta ,lfp_fname,path_name); %(samp0, nSamp, meta, binName, path);

nChan = str2double(LFP.meta.nSavedChans);

nFileSamp = str2double(LFP.meta.fileSizeBytes) / (2 * nChan);
LFP.nSamp = min(LFP.nSamp, nFileSamp - samp0);

sizeA = [nChan, LFP.nSamp];

fid = fopen(lfp_bin_file_path, 'rb');
fseek(fid, samp0 * 2 * nChan, 'bof');
if length(channels) == 384
    % Load all. It's more efficient to do it all at once instead of channel
    % by channel.
    LFP.data_uV = fread(fid, sizeA, 'int16=>single');
    LFP.INFO = TBL;
else
    % load a subset
    LFP.data_uV = nan(length(channels),LFP.nSamp,'single');

    for iCh = 1:length(channels)
        %  Go to start but also skip to the target channel.
        fseek(fid, samp0 * 2 * nChan + channels0(iCh)*2, 'bof');
        % ftell(fid)
        LFP.data_uV(iCh,:) = fread(fid, [1 LFP.nSamp], 'int16=>single',(nChan-1)*2); % remember skip is in bytes to skip, not 'records' or int16s
        fprintf('%d,',iCh)
    end
    LFP.INFO = TBL(channels,:);
end
fclose(fid);

LFP.data_uV = NPXL_Convert_to_uV(LFP.data_uV,LFP.meta);

if PLOT_IT

    t = start_end_rec(1):start_end_rec(2);
    t_min = t/LFP.sFreq/60;
    figure
    subplot(2,1,1)
    imagesc(t_min(1:20:end),[],LFP.data_uV(:,1:20:end));colorbar
    xlabel('min')
    subplot(2,1,2)
    plot(t_min(1:20:end), LFP.data_uV(1:3,1:20:end));
    axis tight
    xlabel('min');ylabel('uV')
end
