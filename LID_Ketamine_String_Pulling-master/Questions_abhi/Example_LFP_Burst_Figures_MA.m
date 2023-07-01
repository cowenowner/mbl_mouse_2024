% examples of LFP filtered trace 80 Hz and ket 50 Hz and some baseline
% examples
LFP_cln = LK_Load_and_Clean_LFP_Abhi(pwd,'M1_Lesion_S6_Ch40.mat');
fq_bands = {'gamma_80' [40 75]}; % or gamma_50

% find a 1 sec interval for 80 Hz or 50 Hz or baseline
ana_range = [-90 30];
times = [-66 -65; -89 -88; -6 -5; -11 -10; 4 5; 8 9];
% t_50 = [];
% t_base = [];
% Make ketamine inj 0
t_sec = (0:(length(LFP_cln.LFP)-1))/LFP_cln.sFreq;
t_uS = LFP_cln.t_uS - LFP.Ket_Start_min*60e6;
t_sec = LFP_cln.t_uS - LFP.Ket_Start_min*1e6;
IX1 = t_uS > ana_range(1)*60e6 & t_uS < ana_range(2)*60e6;
t_uS = t_uS(IX1);
IX2 = t_sec > ana_range(1)*1e6 & t_sec < ana_range(2)*1e6;
t_sec = t_sec(IX2);
F = SPEC_create_filters(fq_bands,LFP_cln.sFreq);
L = filtfilt(F{2},LFP_cln.LFP(IX1));
for ii = 1:Rows(times)
% get the time stamps
% ix1 = binsearch(t_uS, times(ii,1)*60e6);
% ix2 = binsearch(t_uS, times(ii,2)*60e6);
ix1 = binsearch(t_sec, times(ii,1)*1e6);
ix2 = binsearch(t_sec, times(ii,2)*1e6);
% % plotting raw lfp for the given time 
figure
plot(t_sec(ix1:ix2),LFP_cln.LFP(ix1:ix2))
title(num2str(ii))
s = abs(hilbert(L(ix1:ix2)));
%plotting the power of the filtered signal
hold on
plot(t_sec(ix1:ix2),L(ix1:ix2))
% plot(t_uS(ix1:ix2),s)
title(num2str(ii))
end

%% Plotting the burst examples
% Load the relevant TBL from the analysis folder
IX_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control';
IX_HV_LV = IX_CM1 & TBL.LocVar_base > 1.3 & TBL.LocVar_post2 < .9;
Burst = TBL(IX_HV_LV,:);

% D = load('Dset_Rat_315_36.mat');
D = load('Dset_Rat_315_34.mat');
TT = 7;
% Ketamine injection time for Session 36 in 315
ket_inj = LFP.Ket_Start_min; %192 mins
ket_inj = 187;
t_sec = (0:(length(LFP.data)-1))/LFP.LFP_sFreqj;
t_uS = 1e6*t_sec(:);
base_t = [ket_inj-50 ket_inj-40];
ket_t = [ket_inj+3 ket_inj+33];
IX1 = D.SP(TT).t_uS >= base_t(1)*60e6 & D.SP(TT).t_uS <= base_t(2)*60e6;
IX2 = D.SP(TT).t_uS >= ket_t(1)*60e6 & D.SP(TT).t_uS <= ket_t(2)*60e6;
IX_t_b = t_sec >= base_t(1)*1e6 & t_sec <= base_t(2)*1e6;
IX_t_k = t_sec >= ket_t(1)*1e6 & t_sec <= ket_t(2)*1e6;
res_sp_ms = D.SP(TT).t_uS(IX1)/60e6;
[AC,x] = AutoCorrArray(D.SP(TT).t_uS(IX1)/1000,2,150);
figure
plot(x,AC)
title('AutoCorr NRN 7 Session 34 Control M1 Baseline')
[AC,x] = AutoCorrArray(D.SP(TT).t_uS(IX2)/1000,2,150);
figure
plot(x,AC)
title('AutoCorr NRN 7 Session 34 Control M1 Post ket')
figure
subplot(211)
HistISI(D.SP(TT).t_uS(IX1)/100)
title('ISI NRN 7 Session 34 Control M1 Baseline')
subplot(212)
HistISI(D.SP(TT).t_uS(IX2)/100)
title('ISI NRN 7 Session 34 Control M1 Post ket')
figure
histogram(D.SP(TT).t_uS(IX1)/10)

figure;plot((D.SP(TT).t_uS)/1e6,IX2,'.','MarkerSize',80)
ylim([-.5 1.8])
