% Directory that has the LFP subfolder...
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all;
fqs = 1:.2:170;
skip = 1;
win_sec = 40;
data_dir = 'C:\Data\NSandB_course2023\23240\1_7_10_23\vSTR_Ket_Neuro2_g0\vSTR_Ket_Neuro2_g0_imec0';
% data_dir = 'C:\Data\NSandB_course2023\23239\01_7_10_23\vSTR_Ket_Neuro2_g0\vSTR_Ket_Neuro2_g0_imec0';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = dir(fullfile(data_dir,'LFP\*.mat'));

for iF = 1:skip:length(d)
    load(fullfile(data_dir,"LFP",d(iF).name))
    if iF == 1
        DEPTHTBL = NPXL_Depth_From_Meta(LFP.original_meta);
        % figure
        % imagesc(zscore([DEPTHTBL.ChannelID DEPTHTBL.SiteID DEPTHTBL.x_uM DEPTHTBL.y_uM DEPTHTBL.Depth_uM] ))
    end
    % Get the channel and determine the depth.
    ix = find(DEPTHTBL.Ch == LFP.Channel-1);
    Depth_uM = DEPTHTBL.Depth_uM(ix);
    ALL_LFP(iF,:) = single(LFP.Data(:)');
    All_Depth_uM(iF) = Depth_uM;
end
% sort by depth.
[dp,six] = sort(All_Depth_uM);
All_Depth_uM = All_Depth_uM(six);
ALL_LFP = ALL_LFP(six,:);
t_min = (0:(Cols(ALL_LFP)-1))/(LFP.new_sFreq*60);


delta_ALL_LFP = diff(ALL_LFP,[],2);
delta_ALL_LFP = delta_ALL_LFP-mean(delta_ALL_LFP);
delta_delta_ALL_LFP = diff(diff(ALL_LFP,[],2));
delta_delta_ALL_LFP = delta_delta_ALL_LFP-mean(delta_delta_ALL_LFP);
% Pick a deep channel and see if it has delta.
% spectrogram(ALL_LFP(end,:),round(LFP.new_sFreq*win_sec),round(LFP.new_sFreq*win_sec/2),fqs,LFP.new_sFreq,'yaxis');

GIX = t_min > 35 & t_min < 50;
t_min = t_min(GIX);
ALL_LFP = ALL_LFP(:,GIX);
pwelch(ALL_LFP(end-4,:),round(LFP.new_sFreq*win_sec),round(LFP.new_sFreq*win_sec/2),fqs,LFP.new_sFreq);
delta_peak = 3;
theta_peak = 6;
target_ch = Rows(ALL_LFP)-2;
delta_filt = designfilt('bandpassiir','FilterOrder',10, ...
                    'HalfPowerFrequency1',delta_peak-1.5,'HalfPowerFrequency2',delta_peak+1.5, ...
                    'SampleRate',LFP.new_sFreq);

HFO_filt = designfilt('bandpassiir','FilterOrder',16, ...
                    'HalfPowerFrequency1',135,'HalfPowerFrequency2',150, ...
                    'SampleRate',LFP.new_sFreq);
% freqz(HFO_filt)
delta_sig = filtfilt(delta_filt,ALL_LFP(target_ch,:));
delta_env = envelope(abs(hilbert(delta_sig)));
HFO_sig = filtfilt(HFO_filt,ALL_LFP(target_ch,:));
HFO_env = envelope(abs(hilbert(HFO_sig)));
% Find periods of decent delta.
[px,trough_ix] = findpeaks(delta_sig*-1); % invert so that I get troughs.
px = px*-1;
GIX2 = px < -3*std(px);
delta_times_sec = t_min(GIX2)/60;

figure;plot(delta_sig);hold on;plot(trough_ix(GIX2),zeros(size(trough_ix(GIX2))),'r>')

% make a PETH aligned on the trough
% function [M, ix, x_sec] = PETH_EEG_simple(EEG_t_data, alignments_t, samples_before, samples_after,sFreq, PLOT_IT)
for iR = 1:Rows(ALL_LFP)
    HFO_sig = filtfilt(HFO_filt,ALL_LFP(iR,:));
    HFO_env = envelope(abs(hilbert(HFO_sig)));
    [M,ix,x_sec] = PETH_EEG_simple([t_min(:)/60 HFO_env(:)],delta_times_sec,LFP.new_sFreq,LFP.new_sFreq,LFP.new_sFreq);
    MN_PETH_EEG(iR,:) = mean(M);
end
MN_PETH_EEG_norm = Z_scores(MN_PETH_EEG')';
figure
subplot(2,1,1)
imagesc(x_sec,[],MN_PETH_EEG_norm)
subplot(2,1,2)
plot_confidence_intervals(x_sec,MN_PETH_EEG_norm)
plot_vert_line_at_zero

figure
subplot(2,1,1)
imagesc(x_sec,[],diff(MN_PETH_EEG,[],1))
subplot(2,1,2)
plot_confidence_intervals(x_sec,diff(MN_PETH_EEG,[],1))
plot_vert_line_at_zero



figure
imagesc(ALL_LFP(:,:))
figure
imagesc(delta_ALL_LFP(:,:))
figure
imagesc(delta_delta_ALL_LFP(:,:))
caxis([-200,200])

