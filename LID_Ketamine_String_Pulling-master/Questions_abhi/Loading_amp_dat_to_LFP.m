% Parameters specifications
LFP_sFreq = 500;
fs_initial = 30000;
fqs = 1:.5:170;
LFP = [];
t_80 = [0 200];
% Load .dat file
fp = fopen('amp-B-060.dat','rb');
D = fread(fp,'int16');
fclose(fp);
% resample the data to 500 Hz
LFP.data = int16(resample(D,LFP_sFreq, fs_initial));
LFP.to_uV_conversion = 0.195;
LFP.LFP_sFreqj = LFP_sFreq;
LFP.original_sFreq = fs_initial;
LFP.fname = 'amp-B-060.dat';
% Clean the LFP; not sure how much more this helps
LFP_clean =  LK_Load_and_Clean_LFP('E:\LFP_data_Ket_SingleUnit\Rat344\09', 'amp-B-060.mat');
% Get the time stamps for the time of interest
IX = LFP_clean.t_uS > t_80(1)*60e6 & LFP_clean.t_uS < t_80(2)*60e6;

% plotting the data after smoothing 
[S,w,t] = spectrogram(LFP_clean.LFP(IX),LFP_clean.sFreq*5,LFP_clean.sFreq*3,fqs, LFP_clean.sFreq);
t_min = t/60;
S = 10*log10(abs(S))';
%     p = prctile(S,[.5 99.5]);
%     BIX = S<p(1,:);
%     BIX = BIX | S>p(2,:);
%     S(BIX) = nan;
Sn = movmedian(S,30,'omitnan');
Sn = conv_filter(S,hanning(15)/sum(hanning(15))); % smooth data

BASEIX = t_min < 10;
Sb = Sn-nanmean(Sn(BASEIX,:),1);
Sb = Sb./nanstd(Sn(BASEIX,:),[],1);
figure
imagesc(t_min,fqs,Sn');
axis xy
colorbar
caxis(prctile(Sn(:),[1 99]))

% creating 80 Hz filter and plotting filteres signal
F = SPEC_create_filters('gamma_80',LFP_sFreq);
L = filtfilt(F{1},LFP_clean.LFP(IX));
figure
plot(L)

% psd plot
Lp = LFP_clean.LFP(IX);
figure
pwelch(Lp,LFP_clean.sFreq,LFP_clean.sFreq/2,fqs,LFP_clean.sFreq)


