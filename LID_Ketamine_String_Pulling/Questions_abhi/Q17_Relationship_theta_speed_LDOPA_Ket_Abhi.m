function OUT = Q17_Relationship_theta_speed_LDOPA_Ket_Abhi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abhi 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OUT = [];

[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();

PLOT_IT = false;
% binsize_Q_S_ms = 500; % IMU Speed
% bin_size_ac_ms = 2;
% acorr_width_ms = 200;

% ket_gamma_fqs = 35:.5:75;
theta_filt = [5 11];
spec_fq = [1:.25:100];

speed_bin_cm_sec = 1;

% OUT.binsize_Q_S_ms = binsize_Q_S_ms;

OUT.SP = SP;
OUT.DEPTHS = DEPTHS; % Send this info out of the function for meta analysis.
OUT.META = META;
OUT.EVT = EVT;

OUT.theta_filt = theta_filt;
OUT.speed_bin_cm_sec = speed_bin_cm_sec;
OUT.spec_fq = spec_fq;

e_start_end_uS = [];

if any(E.EventID == 'KetInjectionStart') 
        e_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
        e_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-55 -35; -25 -5; 2 12; 17 27; 60 90];
        big_peri_event_min = [-55 30];

elseif any(E.EventID == 'LDOPAInjectionStart')
        e_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
        e_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 25 55; 65 75; 80 90; 120 140];
        big_peri_event_min = [-25 90];
end

% return
OUT.intervals_around_evt_min = intervals_around_evt_min;
TIMES.EventStartEndUsec(1,1) = e_start_end_uS(1);
TIMES.EventStartEndUsec(1,2) = e_start_end_uS(2);
TIMES.BaselineUsec(1,1) = TIMES.EventStartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.BaselineUsec(1,2) = TIMES.EventStartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.PostInjectionUsec(1,1) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(2,1)*60*1e6;
TIMES.PostInjectionUsec(1,2) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(2,2)*60*1e6;
TIMES.PostInjectionUsec(2,1) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(3,1)*60*1e6;
TIMES.PostInjectionUsec(2,2) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(3,2)*60*1e6;
TIMES.PostInjectionUsec(3,1) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(4,1)*60*1e6;
TIMES.PostInjectionUsec(3,2) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(4,2)*60*1e6;
TIMES.PostInjectionUsec(4,1) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(5,1)*60*1e6;
TIMES.PostInjectionUsec(4,2) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(5,2)*60*1e6;


TIMES.PeriKetamineUsec = [TIMES.EventStartEndUsec(1) + big_peri_event_min(1)*60*1e6 TIMES.EventStartEndUsec(1) + big_peri_event_min(2)*60*1e6 ];

OUT.TM = TIMES;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load inertial, POS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IMU = LK_Load_and_Process_IMU('Inertial_data.mat');
POS = LK_Load_and_Clean_POS_Abhi('POS.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find indices for the times for baseline, and the two post-ketamine periods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IX1_P = POS.Time_uS > TIMES.BaselineUsec(1) & POS.Time_uS < TIMES.BaselineUsec(2);
IX2_P = POS.Time_uS > TIMES.PostInjectionUsec(1,1) & POS.Time_uS < TIMES.PostInjectionUsec(1,2);
IX3_P = POS.Time_uS > TIMES.PostInjectionUsec(2,1) & POS.Time_uS < TIMES.PostInjectionUsec(2,2);
IX4_P = POS.Time_uS > TIMES.PostInjectionUsec(3,1) & POS.Time_uS < TIMES.PostInjectionUsec(3,2);
% IX5_P = POS.Time_uS > TIMES.PeriKetamineUsec(1) & POS.Time_uS < TIMES.PeriKetamineUsec(2);
% OUT.Speed_Base_Sal_Early_Late = [nanmean(POS.speed(IX1_P)) nanmean(POS.speed(IX2_P)) nanmean(POS.speed(IX3_P)) nanmean(POS.speed(IX4_P))];

IX1 = IMU.t_uS > TIMES.BaselineUsec(1) & IMU.t_uS < TIMES.BaselineUsec(2);
IX2 = IMU.t_uS > TIMES.PostInjectionUsec(1,1) & IMU.t_uS < TIMES.PostInjectionUsec(1,2);
IX3 = IMU.t_uS > TIMES.PostInjectionUsec(2,1) & IMU.t_uS < TIMES.PostInjectionUsec(2,2);
IX4 = IMU.t_uS > TIMES.PostInjectionUsec(3,1) & IMU.t_uS < TIMES.PostInjectionUsec(3,2);
% IX5 = IMU.t_uS > TIMES.PeriKetamineUsec(1) & IMU.t_uS < TIMES.PeriKetamineUsec(2);
% Store information for analysis
% OUT.Jerk_base_post1_post2_post3 = [nanmean(IMU.absjerk(IX1)) nanmean(IMU.absjerk(IX2)) nanmean(IMU.absjerk(IX3)) nanmean(IMU.absjerk(IX4))];
% OUT.JerkPC1_base_post1_post2_post3 = [nanmean(IMU.absjerkpc(IX1)) nanmean(IMU.absjerkpc(IX2)) nanmean(IMU.absjerkpc(IX3)) nanmean(IMU.absjerkpc(IX4))];
% OUT.IMU_Speed_base_post1_post2_post3 = [nanmean(IMU.speed(IX1)) nanmean(IMU.speed(IX2)) nanmean(IMU.speed(IX3)) nanmean(IMU.speed(IX4))];

IMU_Speed_base = [IMU.t_uS(IX1) IMU.speed(IX1)];
IMU_Speed_post1 = [IMU.t_uS(IX2) IMU.speed(IX2)];
IMU_Speed_post2 = [IMU.t_uS(IX3) IMU.speed(IX3)];
IMU_Speed_post3 = [IMU.t_uS(IX4) IMU.speed(IX4)];
% IMU_Speed_All = [IMU.t_uS(IX5) IMU.speed(IX5)];

POS_Speed_base = [POS.Time_uS(IX1_P) POS.speed(IX1_P)];
POS_Speed_post1 = [POS.Time_uS(IX2_P) POS.speed(IX2_P)];
POS_Speed_post2 = [POS.Time_uS(IX3_P) POS.speed(IX3_P)];
POS_Speed_post3 = [POS.Time_uS(IX4_P) POS.speed(IX4_P)];

OUT.IMU_Speed_base = IMU_Speed_base(:,2);
OUT.IMU_Speed_base_t_uS = IMU_Speed_base(:,1);
OUT.IMU_Speed_post1 = IMU_Speed_post1(:,2);
OUT.IMU_Speed_post1_t_uS = IMU_Speed_post1(:,1);
OUT.IMU_Speed_post2 = IMU_Speed_post2(:,2);
OUT.IMU_Speed_post2_t_uS = IMU_Speed_post2(:,1);
OUT.IMU_Speed_post3 = IMU_Speed_post3(:,2);
OUT.IMU_Speed_post3_t_uS = IMU_Speed_post3(:,1);

OUT.POS_Speed_base = POS_Speed_base(:,2);
OUT.POS_Speed_base_t_uS = POS_Speed_base(:,1);
OUT.POS_Speed_post1 = POS_Speed_post1(:,2);
OUT.POS_Speed_post1_t_uS = POS_Speed_post1(:,1);
OUT.POS_Speed_post2 = POS_Speed_post2(:,2);
OUT.POS_Speed_post2_t_uS = POS_Speed_post2(:,1);
OUT.POS_Speed_post3 = POS_Speed_post3(:,2);
OUT.POS_Speed_post3_t_uS = POS_Speed_post3(:,1);

% POS_Speed_All = [POS.Time_uS(IX5) POS.speed(IX5)];
% 
% % Bin the POS and IMU speed data
% % Create start and end times
% edges_b = TIMES.BaselineUsec(1,1):binsize_Q_S_ms*1000:TIMES.BaselineUsec(1,2);
% [IMUbase,bins_uS] = Bin_IMU_data(IMU_Speed_base, edges_b);
% edges_p1 = TIMES.PostInjectionUsec(1,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(1,2);
% [IMUpost1,y1] = Bin_IMU_data(IMU_Speed_post1, edges_p1);
% edges_p2 = TIMES.PostInjectionUsec(2,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(2,2);
% [IMUpost2,y2] = Bin_IMU_data(IMU_Speed_post2, edges_p2);
% edges_p3 = TIMES.PostInjectionUsec(3,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(3,2);
% [IMUpost3,y3] = Bin_IMU_data(IMU_Speed_post3, edges_p3);
% 
% edges_all = TIMES.PeriKetamineUsec(1):binsize_Q_S_ms*1000:TIMES.PeriKetamineUsec(2);
% [IMUall,IMUall_x_uS] = Bin_IMU_data(IMU_Speed_All, edges_all);
% 
% OUT.IMU_Speed_base = IMUbase;
% OUT.IMU_base_bins_uS = bins_uS(:,1) + binsize_Q_S_ms*1000/2;
% OUT.IMU_Speed_post1 = IMUpost1;
% OUT.IMU_psot1_bins_uS = y1(:,1) + binsize_Q_S_ms*1000/2;
% OUT.IMU_Speed_post2 = IMUpost2;
% OUT.IMU_post2_bins_uS = y2(:,1) + binsize_Q_S_ms*1000/2;
% OUT.IMU_Speed_post3 = IMUpost3;
% OUT.IMU_post3_bins_uS = y3(:,1) + binsize_Q_S_ms*1000/2;
% OUT.IMU_Speed_All = IMUall;
% OUT.IMU_All_bins_uS = IMUall_x_uS(:,1) + binsize_Q_S_ms*1000/2;
% % 
% % IMU_post2_bins_uS = OUT.IMU_post2_bins_uS;
% % IMU_post2_bins_sec = IMU_post2_bins_uS/1e6;
% % IMU_post2_sec = IMU_post2_bins_sec - IMU_post2_bins_sec(1);
% 
% [POSbase,bins_uS_p] = Bin_IMU_data(POS_Speed_base, edges_b);
% [POSpost1,y1_p] = Bin_IMU_data(POS_Speed_post1, edges_p1);
% [POSpost2,y2_p] = Bin_IMU_data(POS_Speed_post2, edges_p2);
% [POSpost3,y3_p] = Bin_IMU_data(POS_Speed_post3, edges_p3);
% % [POSall,POSall_x_uS] = Bin_IMU_data(POS_Speed_All, edges_all);
% 
% OUT.POS_Speed_base = POSbase;
% OUT.POS_base_bins_uS = bins_uS_p(:,1) + binsize_Q_S_ms*1000/2;
% OUT.POS_Speed_post1 = POSpost1;
% OUT.POS_psot1_bins_uS = y1_p(:,1) + binsize_Q_S_ms*1000/2;
% OUT.POS_Speed_post2 = POSpost2;
% OUT.POS_post2_bins_uS = y2_p(:,1) + binsize_Q_S_ms*1000/2;
% OUT.POS_Speed_post3 = POSpost3;
% OUT.POS_post3_bins_uS = y3_p(:,1) + binsize_Q_S_ms*1000/2;
% OUT.POS_Speed_All = POSall;
% OUT.POS_All_bins_uS = POSall_x_uS(:,1) + binsize_Q_S_ms*1000/2;
% 
% POS_post2_bins_uS = OUT.POS_post2_bins_uS;
% POS_post2_bins_sec = POS_post2_bins_uS/1e6;
% POS_post2_sec = POS_post2_bins_sec - POS_post2_bins_sec(1);


LFPfiles = find_files(fullfile('M1*.mat'));
LFP = LK_Load_and_Clean_LFP_Abhi('',LFPfiles{1});
L = [];
Base_IX = LFP.t_uS > TIMES.BaselineUsec(1) & LFP.t_uS < TIMES.BaselineUsec(2);
Ket0_IX = LFP.t_uS > TIMES.PostInjectionUsec(1,1) & LFP.t_uS < TIMES.PostInjectionUsec(1,2);
Ket1_IX = LFP.t_uS > TIMES.PostInjectionUsec(2,1) & LFP.t_uS < TIMES.PostInjectionUsec(2,2);
Ket2_IX = LFP.t_uS > TIMES.PostInjectionUsec(3,1) & LFP.t_uS < TIMES.PostInjectionUsec(3,2);
Ket3_IX = LFP.t_uS > TIMES.PostInjectionUsec(4,1) & LFP.t_uS < TIMES.PostInjectionUsec(4,2);
% get psd 
psdb = pwelch(LFP.LFP(Base_IX),LFP.sFreq,LFP.sFreq/2,spec_fq,LFP.sFreq);
psdb = 10*log10(abs(psdb));
psd0 = pwelch(LFP.LFP(Ket0_IX),LFP.sFreq,LFP.sFreq/2,spec_fq,LFP.sFreq);
psd0 = 10*log10(abs(psd0));
psd1 = pwelch(LFP.LFP(Ket1_IX),LFP.sFreq,LFP.sFreq/2,spec_fq,LFP.sFreq);
psd1 = 10*log10(abs(psd1));
psd2 = pwelch(LFP.LFP(Ket2_IX),LFP.sFreq,LFP.sFreq/2,spec_fq,LFP.sFreq);
psd2 = 10*log10(abs(psd2));
psd3 = pwelch(LFP.LFP(Ket3_IX),LFP.sFreq,LFP.sFreq/2,spec_fq,LFP.sFreq);
psd3 = 10*log10(abs(psd3));

OUT.PSD_base = psdb;
OUT.PSD_ket0 = psd0;
OUT.PSD_ket1 = psd1;
OUT.PSD_ket2 = psd2;
OUT.PSD_ket3 = psd3;

F = SPEC_create_filters(theta_filt,LFP.sFreq);
L.LFP_ket0_t_uS = LFP.t_uS(Ket0_IX);
L.LFP_ket0_filt = filtfilt(F{1},LFP.LFP(Ket0_IX));
L.LFP_ket1_t_uS = LFP.t_uS(Ket1_IX);
L.LFP_ket1_filt = filtfilt(F{1},LFP.LFP(Ket1_IX));
L.LFP_ket2_t_uS = LFP.t_uS(Ket2_IX);
L.LFP_ket2_filt = filtfilt(F{1},LFP.LFP(Ket2_IX));    

OUT.LFP_sfreq = LFP.sFreq;

%%
% [S,f,t] = spectrogram(L.LFP_ket1_filt,LFP.sFreq*spec_win_sec,LFP.sFreq*(spec_win_sec/2),ket_gamma_fqs, LFP.sFreq);
% sFreq_of_spectrogram = 1/median(diff(t));
% t_min = t/60;
% S = 10*log10(abs(S))';
% 
% han_size_sec = 40;
% han_size_points = han_size_sec*sFreq_of_spectrogram;
% Sn = conv_filter(S,hanning(han_size_points)/sum(hanning(han_size_points)));
%     
% S_norm = Sn-nanmean(Sn,1);
% S_norm = S_norm./nanstd(Sn,[],1);
% 
% OUT.S_norm = S_norm;
% OUT.ket_gamma_fqs = f;
% OUT.S_t_sec = t;
% 
% if PLOT_IT
%     figure
%     imagesc(t_min,ket_gamma_fqs,S_norm');
%     
%     axis xy
%     colorbar
%     caxis(prctile(S_norm(:),[1 98.6]))
%     xlabel('Time min')
%     ylabel('HZ')
%     title('Spectrogram post ket 2-12 min Rat 352 Ses 5')
%     
%     figure
%     plot(POSpost2)
% end
% 
% Spec_interp = interp1(t, S_norm, POS_post2_sec, 'makima');
% Spec_interp_IMU = interp1(t, S_norm, IMU_post2_sec, 'makima');
% 
% OUT.Spec_interp = Spec_interp;
% 
% if PLOT_IT
%     figure
%     imagesc(POSpost2, ket_gamma_fqs, Spec_interp')
%     axis xy
%     colorbar
%     caxis(prctile(Spec_interp(:),[1 98.6]))
%     
%     figure
%     imagesc(POSpost2, ket_gamma_fqs, Spec_interp_IMU')
%     axis xy
%     colorbar
%     caxis(prctile(Spec_interp(:),[1 98.6]))
% end
% 
% [B, I] = sort(POSpost2);
% wtf = sort_matrix(Spec_interp, 'by_value', I);
% 
% [B_IMU, I_IMU] = sort(IMUpost2);
% wtf_IMU = sort_matrix(Spec_interp_IMU, 'by_value', I_IMU);
% % wtf = sort_matrix(Spec_interp,'by_value',POSpost2);
% if PLOT_IT
% 
%     figure; plot(wtf)
%     figure; imagesc(B, ket_gamma_fqs, wtf')
%     axis xy
%     colorbar
%     xlabel('POS speed cm/s')
%     ylabel('HZ')
%     title('Speed vs frequency post ket 2-30 min Rat 352 Ses 5')
%     
%     figure; imagesc(B_IMU, ket_gamma_fqs, wtf_IMU')
%     axis xy
%     colorbar
%     xlabel('IMU speed')
%     ylabel('HZ')
%     title('Speed vs frequency post ket 2-30 min Rat 352 Ses 5')
% end

%% Before ket either peak LDOPA or Post saline
% Get the power
p0 = abs(hilbert(L.LFP_ket0_filt));
p0(p0>90) = nan;
Power_gamma_ket0 = movmean(p0,LFP.sFreq/2,'omitnan');
Power_interp_ket0 = interp1(L.LFP_ket0_t_uS, Power_gamma_ket0, POS_Speed_post1(:,1));

% Instantaneous freq
v0 = instfreq(L.LFP_ket0_filt,LFP.sFreq,'Method','hilbert');

v0(end+1) = v0(end);

v0(abs(v0)>30) = nan;

Inst_freq_gamma_ket0 = movmean(v0,LFP.sFreq/2,'omitnan');

Inst_freq_interp_ket0 = interp1(L.LFP_ket0_t_uS, Inst_freq_gamma_ket0, POS_Speed_post1(:,1));

% get start and end times of the 2 cm/s bins for POS speed
Speed_bins_ket0 = min(POS_Speed_post1(:,2)):speed_bin_cm_sec:max(POS_Speed_post1(:,2));
min_bin_samples_speed = 100;
Bin_freq_by_speed_ket0 = [];
Bin_power_by_speed_ket0 = [];
for ii = 1:length(Speed_bins_ket0)-1
    IX_speed_bin = POS_Speed_post1(:,2) >= Speed_bins_ket0(ii) & POS_Speed_post1(:,2) <= Speed_bins_ket0(ii+1);
    if sum(IX_speed_bin) > min_bin_samples_speed
        Bin_freq_by_speed_ket0(ii) = mean(Inst_freq_interp_ket0(IX_speed_bin));
        Bin_power_by_speed_ket0(ii) = mean(Power_interp_ket0(IX_speed_bin));
    else
        Bin_freq_by_speed_ket0(ii) = nan;
        Bin_power_by_speed_ket0(ii) = nan;
    end
end
Speed_bin_center_ket0 = Speed_bins_ket0 + speed_bin_cm_sec/2;
Speed_bin_center_ket0 = Speed_bin_center_ket0(1:end-1);

OUT.Bin_power_by_speed_ket0 = Bin_power_by_speed_ket0;
OUT.Bin_freq_by_speed_ket0 = Bin_freq_by_speed_ket0;
OUT.Speed_bin_center_ket0 = Speed_bin_center_ket0;

% Ketamine first half
% Get the power
p1 = abs(hilbert(L.LFP_ket1_filt));
p1(p1>90) = nan;
Power_gamma_ket1 = movmean(p1,LFP.sFreq/2,'omitnan');
Power_interp_ket1 = interp1(L.LFP_ket1_t_uS, Power_gamma_ket1, POS_Speed_post2(:,1));

% Instantaneous freq
v1 = instfreq(L.LFP_ket1_filt,LFP.sFreq,'Method','hilbert');

v1(end+1) = v1(end);

v1(abs(v1)>30) = nan;

Inst_freq_gamma_ket1 = movmean(v1,LFP.sFreq/2,'omitnan');

Inst_freq_interp_ket1 = interp1(L.LFP_ket1_t_uS, Inst_freq_gamma_ket1, POS_Speed_post2(:,1));

% get start and end times of the 2 cm/s bins for POS speed
Speed_bins_ket1 = min(POS_Speed_post2(:,2)):speed_bin_cm_sec:max(POS_Speed_post2(:,2));
min_bin_samples_speed = 100;
Bin_freq_by_speed_ket1 = [];
Bin_power_by_speed_ket1 = [];
for ii = 1:length(Speed_bins_ket1)-1
    IX_speed_bin = POS_Speed_post2(:,2) >= Speed_bins_ket1(ii) & POS_Speed_post2(:,2) <= Speed_bins_ket1(ii+1);
    if sum(IX_speed_bin) > min_bin_samples_speed
        Bin_freq_by_speed_ket1(ii) = mean(Inst_freq_interp_ket1(IX_speed_bin));
        Bin_power_by_speed_ket1(ii) = mean(Power_interp_ket1(IX_speed_bin));
    else
        Bin_freq_by_speed_ket1(ii) = nan;
        Bin_power_by_speed_ket1(ii) = nan;
    end
end
Speed_bin_center_ket1 = Speed_bins_ket1 + speed_bin_cm_sec/2;
Speed_bin_center_ket1 = Speed_bin_center_ket1(1:end-1);

OUT.Bin_power_by_speed_ket1 = Bin_power_by_speed_ket1;
OUT.Bin_freq_by_speed_ket1 = Bin_freq_by_speed_ket1;
OUT.Speed_bin_center_ket1 = Speed_bin_center_ket1;

if PLOT_IT
    figure
    plot(Speed_bin_center_ket1, Bin_freq_by_speed_ket1)
end

%% Ketamine second half 

p2 = abs(hilbert(L.LFP_ket2_filt));
p2(p2>90) = nan;
Power_gamma_ket2 = movmean(p2,LFP.sFreq/2,'omitnan');
Power_interp_ket2 = interp1(L.LFP_ket2_t_uS, Power_gamma_ket2, POS_Speed_post3(:,1));

v2 = instfreq(L.LFP_ket2_filt,LFP.sFreq,'Method','hilbert');

v2(end+1) = v2(end);

v2(abs(v2)>30) = nan;

Inst_freq_gamma_ket2 = movmean(v2,LFP.sFreq/2,'omitnan');

Inst_freq_interp_ket2 = interp1(L.LFP_ket2_t_uS, Inst_freq_gamma_ket2, POS_Speed_post3(:,1));

% get start and end times of the 2 cm/s bins for POS speed
Speed_bins_ket2 = min(POS_Speed_post3(:,2)):speed_bin_cm_sec:max(POS_Speed_post3(:,2));
Bin_freq_by_speed_ket2 = [];
Bin_power_by_speed_ket2 = [];

for ii = 1:length(Speed_bins_ket2)-1
    IX_speed_bin = POS_Speed_post3(:,2) >= Speed_bins_ket2(ii) & POS_Speed_post3(:,2) <= Speed_bins_ket2(ii+1);
    if sum(IX_speed_bin) > min_bin_samples_speed
        Bin_freq_by_speed_ket2(ii) = mean(Inst_freq_interp_ket2(IX_speed_bin));
        Bin_power_by_speed_ket2(ii) = mean(Power_interp_ket2(IX_speed_bin));
    else
        Bin_freq_by_speed_ket2(ii) = nan;
        Bin_power_by_speed_ket2(ii) = nan;
    end
end
Speed_bin_center_ket2 = Speed_bins_ket2 + speed_bin_cm_sec/2;
Speed_bin_center_ket2 = Speed_bin_center_ket2(1:end-1);

OUT.Bin_power_by_speed_ket2 = Bin_power_by_speed_ket2;
OUT.Bin_freq_by_speed_ket2 = Bin_freq_by_speed_ket2;
OUT.Speed_bin_center_ket2 = Speed_bin_center_ket2;
OUT.min_bin_samples_speed = min_bin_samples_speed;