function [NIDQ] = NPXL_Extract_NIDQ(nidq_bin_file_path, channels, start_end_rec )
% function [NIDQ] = NPXL_Extract_LFP(nidq_bin_file_path, channels, start_end_rec )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nidq_bin_file_path = name of the .nidq.
% channels = Channels are in 1 based numbering so channel 0 on Imec = channel 1 here.
% start_end_rec = start_end_rec are also 1 based (so first record is 1, not 0).
%
% Outputs NIDQ structure 
%
% Cowen 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
PLOT_IT = false;
if nargin <2
    channels = 1:7;
end

[path_name, tmp] = fileparts(nidq_bin_file_path);
lfp_fname = [tmp '.bin'];
meta_fname = strrep(lfp_fname,'.nidq.bin','.nidq.meta');

obj = SGLX_Class;
NIDQ.meta = obj.ReadMeta(meta_fname,path_name);
NIDQ.nChan = str2double(NIDQ.meta.nSavedChans);
NIDQ.sFreq = str2double(NIDQ.meta.niSampRate);
NIDQ.channels = channels;
NIDQ.nFileSamp = str2double(NIDQ.meta.fileSizeBytes) / (2 * NIDQ.nChan);
NIDQ.duration_of_recording_sec = NIDQ.nFileSamp/NIDQ.sFreq;
NIDQ.fI2V = str2double(NIDQ.meta.niAiRangeMax) / 32768;
% dataVolts = dataInt * fI2V / gain.

if nargin <3
    start_end_rec = [0 NIDQ.nFileSamp-NIDQ.nChan*2];
end

samp0 = start_end_rec(1) -1;
samp0 = max(samp0, 0);

NIDQ.nSamp = start_end_rec(2) - start_end_rec(1) + 1;
NIDQ.start_rec = start_end_rec(1);

nFileSamp = str2double(NIDQ.meta.fileSizeBytes) / (2 * NIDQ.nChan);
NIDQ.nSamp = min(NIDQ.nSamp, nFileSamp - samp0);

sizeA = [NIDQ.nChan, NIDQ.nSamp];

fid = fopen(nidq_bin_file_path, 'rb');
fseek(fid, samp0 * 2 * NIDQ.nChan, 'bof');
NIDQ.data_V = fread(fid, sizeA, 'int16=>single');
fclose(fid);

NIDQ.data_V = obj.GainCorrectNI(NIDQ.data_V, 1:7, NIDQ.meta);

t = start_end_rec(1):start_end_rec(2);
NIDQ.t_sec = t/NIDQ.sFreq;

NIDQ.data_V = NIDQ.data_V(channels,:);

if PLOT_IT
    figure
    subplot(2,1,1)
    imagesc(NIDQ.t_sec(1:20:end)/60,[],NIDQ.data_V(:,1:20:end));colorbar
    xlabel('min')
    subplot(2,1,2)
    plot  (NIDQ.t_sec(1:20:end)/60, NIDQ.data_V(:,1:20:end));
    axis tight
    xlabel('min');ylabel('V')
end
