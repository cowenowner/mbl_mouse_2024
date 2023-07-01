function OUT = Q13_Does_FR_correlate_with_movement_during_ketamine_Abhi()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abhi 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OUT = [];

binsize_Q_S_ms = 500; % 0.5 sec for FR and IMU Speed
bin_size_ac_ms = 2;
acorr_width_ms = 200;
PLOT_IT = false;
% Post_Ket_time = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function loads a ton of things that are typically used
% for any analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;

OUT.binsize_Q_S_ms = binsize_Q_S_ms;

OUT.SP = SP;
OUT.DEPTHS = DEPTHS; % Send this info out of the function for meta analysis.
OUT.META = META;
OUT.EVT = EVT;

e_start_end_uS = [];

if any(E.EventID == 'KetInjectionStart') && any(E.EventID == 'SalineInjectionStart')
        e_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
        e_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-55 -35; -25 -5; 2 30; 60 80];
        big_peri_event_min = [-55 80];

elseif any(E.EventID == 'LDOPAInjectionStart')
        e_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
        e_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 25 55; 65 90; 120 140];
        big_peri_event_min = [-25 140];
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
IX1 = POS.Time_uS > TIMES.BaselineUsec(1) & POS.Time_uS < TIMES.BaselineUsec(2);
IX2 = POS.Time_uS > TIMES.PostInjectionUsec(1,1) & POS.Time_uS < TIMES.PostInjectionUsec(1,2);
IX3 = POS.Time_uS > TIMES.PostInjectionUsec(2,1) & POS.Time_uS < TIMES.PostInjectionUsec(2,2);
IX4 = POS.Time_uS > TIMES.PostInjectionUsec(3,1) & POS.Time_uS < TIMES.PostInjectionUsec(3,2);
OUT.Speed_Base_Sal_Early_Late = [nanmean(POS.speed(IX1)) nanmean(POS.speed(IX2)) nanmean(POS.speed(IX3)) nanmean(POS.speed(IX4))];

IX1 = IMU.t_uS > TIMES.BaselineUsec(1) & IMU.t_uS < TIMES.BaselineUsec(2);
IX2 = IMU.t_uS > TIMES.PostInjectionUsec(1,1) & IMU.t_uS < TIMES.PostInjectionUsec(1,2);
IX3 = IMU.t_uS > TIMES.PostInjectionUsec(2,1) & IMU.t_uS < TIMES.PostInjectionUsec(2,2);
IX4 = IMU.t_uS > TIMES.PostInjectionUsec(3,1) & IMU.t_uS < TIMES.PostInjectionUsec(3,2);
IX5 = IMU.t_uS > TIMES.PeriKetamineUsec(1) & IMU.t_uS < TIMES.PeriKetamineUsec(2);
% Store information for analysis
% OUT.Jerk_base_post1_post2_post3 = [nanmean(IMU.absjerk(IX1)) nanmean(IMU.absjerk(IX2)) nanmean(IMU.absjerk(IX3)) nanmean(IMU.absjerk(IX4))];
% OUT.JerkPC1_base_post1_post2_post3 = [nanmean(IMU.absjerkpc(IX1)) nanmean(IMU.absjerkpc(IX2)) nanmean(IMU.absjerkpc(IX3)) nanmean(IMU.absjerkpc(IX4))];
OUT.IMU_Speed_base_post1_post2_post3 = [nanmean(IMU.speed(IX1)) nanmean(IMU.speed(IX2)) nanmean(IMU.speed(IX3)) nanmean(IMU.speed(IX4))];

IMU_Speed_base = [IMU.t_uS(IX1) IMU.speed(IX1)];
IMU_Speed_post1 = [IMU.t_uS(IX2) IMU.speed(IX2)];
IMU_Speed_post2 = [IMU.t_uS(IX3) IMU.speed(IX3)];
IMU_Speed_post3 = [IMU.t_uS(IX4) IMU.speed(IX4)];
IMU_Speed_All = [IMU.t_uS(IX5) IMU.speed(IX5)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Always visualize.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movement
if PLOT_IT
    figure
    plot(IMU.t_uS/60e6,IMU.absjerk)
    hold on
    plot(IMU.t_uS/60e6,IMU.speed,'k')
    
    ylabel('Jerk')
    xlabel('min')
    plot_markers_simple(TIMES.EventStartEndUsec/60e6,[],[],GP.Colors.Ketamine)
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'g')
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'c')
    
    yyaxis right
    plot(POS.Time_uS/60e6,POS.speed)
    ylabel('speed from x and y')
    pubify_figure_axis
    axis tight
    title('Position Data')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % speed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    subplot(2,1,1)
    bar(OUT.Speed_Base_Sal_Early_Late)
    set(gca,'XTickLabel',{'before' 'after 1' 'after 2' 'after 3'})
    ylabel('Speed')
    pubify_figure_axis

    subplot(2,1,2)
    bar(OUT.IMU_Speed_base_post1_post2_post3)
    set(gca,'XTickLabel',{'before' 'after 1' 'after 2' 'after 3'})
    ylabel('IMU Speed')
    pubify_figure_axis
end

% Bin the IMU speed data
% Create start and end times
edges = TIMES.BaselineUsec(1,1):binsize_Q_S_ms*1000:TIMES.BaselineUsec(1,2);
[IMUbase,bins_uS] = Bin_IMU_data(IMU_Speed_base, edges);
edges = TIMES.PostInjectionUsec(1,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(1,2);
[IMUpost1,y1] = Bin_IMU_data(IMU_Speed_post1, edges);
edges = TIMES.PostInjectionUsec(2,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(2,2);
[IMUpost2,y2] = Bin_IMU_data(IMU_Speed_post2, edges);
edges = TIMES.PostInjectionUsec(3,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(3,2);
[IMUpost3,y3] = Bin_IMU_data(IMU_Speed_post3, edges);

edges = TIMES.PeriKetamineUsec(1):binsize_Q_S_ms*1000:TIMES.PeriKetamineUsec(2);
[IMUall,IMUall_x_uS] = Bin_IMU_data(IMU_Speed_All, edges);

OUT.IMU_Speed_base = IMUbase;
OUT.IMU_base_bins_uS = bins_uS(:,1) + binsize_Q_S_ms*1000/2;
OUT.IMU_Speed_post1 = IMUpost1;
OUT.IMU_psot1_bins_uS = y1(:,1) + binsize_Q_S_ms*1000/2;
OUT.IMU_Speed_post2 = IMUpost2;
OUT.IMU_post2_bins_uS = y2(:,1) + binsize_Q_S_ms*1000/2;
OUT.IMU_Speed_post3 = IMUpost3;
OUT.IMU_post3_bins_uS = y3(:,1) + binsize_Q_S_ms*1000/2;
OUT.IMU_Speed_All = IMUall;
OUT.IMU_All_bins_uS = IMUall_x_uS(:,1) + binsize_Q_S_ms*1000/2;


% plot it
if PLOT_IT
    figure
    subplot(411)
    plot(IMUbase)
    subplot(412)
    plot(IMUpost1)
    subplot(413)
    plot(IMUpost2)
    subplot(414)
    plot(IMUpost3)
end

%% Single unit activity binning
% binsize_Q_S_ms = 0;
edges = TIMES.BaselineUsec(1,1):binsize_Q_S_ms*1000:TIMES.BaselineUsec(1,2);
[Qbase,bins_uS] = Bin_ts_array(TS, edges);
edges = TIMES.PostInjectionUsec(1,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(1,2);
[Qpost1,x1] = Bin_ts_array(TS, edges);
edges = TIMES.PostInjectionUsec(2,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(2,2);
[Qpost2,x2] = Bin_ts_array(TS, edges);
edges = TIMES.PostInjectionUsec(3,1):binsize_Q_S_ms*1000:TIMES.PostInjectionUsec(3,2);
[Qpost3,x3] = Bin_ts_array(TS, edges);

edges = TIMES.PeriKetamineUsec(1):binsize_Q_S_ms*1000:TIMES.PeriKetamineUsec(2);
[Qall,Qall_x_uS] = Bin_ts_array(TS, edges);

TSaligned = cell(size(TS));
for ii = 1:length(TS)
    if ~isempty(TS{ii})
        TSaligned{ii} = TS{ii} - TIMES.EventStartEndUsec(1);
    end
end
OUT.TSaligned = TSaligned;

bins_s = binsize_Q_S_ms/1000;
mbase = mean(Qbase/bins_s);
m1 = mean(Qpost1/bins_s);
m2 = mean(Qpost2/bins_s);
m3 = mean(Qpost3/bins_s);

OUT.FRates_base = mbase;
OUT.FRates_post1 = m1;
OUT.FRates_post2 = m2;
OUT.FRates_post3 = m3;

OUT.Qbase = Qbase;
OUT.Qbase_bins_uS = bins_uS(:,1) + binsize_Q_S_ms*1000/2;
OUT.Qpost1 = Qpost1;
OUT.Qpost1_bins_uS = x1(:,1) + binsize_Q_S_ms*1000/2;
OUT.Qpost2 = Qpost2;
OUT.Qpost2_bins_uS = x2(:,1) + binsize_Q_S_ms*1000/2;
OUT.Qpost3 = Qpost3;
OUT.Qpost3_bins_uS = x3(:,1) + binsize_Q_S_ms*1000/2;
OUT.Qall = Qall;
OUT.Qall_x_uS = Qall_x_uS(:,1) + binsize_Q_S_ms*1000/2;

% plot and look at data
[Q,start_end_times] = Bin_ts_array(TS, binsize_Q_S_ms*1000); % workhorse function for neuron x time matrix
IX = start_end_times(:,1)<TIMES.EventStartEndUsec(1);
Qn = Q-mean(Q(IX,:));
[~,sc] = pca(Qn);

% Creating groups for plotting in gscatter for PCA 
test = Qall_x_uS(:,1) > TIMES.BaselineUsec(1) & Qall_x_uS(:,1) <= TIMES.BaselineUsec(2);
test1 = Qall_x_uS(:,1) > TIMES.PostInjectionUsec(1,1) & Qall_x_uS(:,1) <= TIMES.PostInjectionUsec(1,2);
test2 = Qall_x_uS(:,1) > TIMES.PostInjectionUsec(2,1) & Qall_x_uS(:,1) <= TIMES.PostInjectionUsec(2,2);
test3 = Qall_x_uS(:,1) > TIMES.PostInjectionUsec(3,1) & Qall_x_uS(:,1) <= TIMES.PostInjectionUsec(3,2);

test = double(test);
test(test1) = 2;
test(test2) = 3;
test(test3) = 4;
drug_g = test;

OUT.Qall_drug_times = drug_g;

% Correlation between IMU speed and neuron activity
[rbase, pbase] = corr(IMUbase,Qbase);
OUT.Corr_speed_FR_r_p_base = [rbase; pbase];
[rpost1, ppost1] = corr(IMUpost1,Qpost1);
OUT.Corr_speed_FR_r_p_post1 = [rpost1; ppost1];
[rpost2, ppost2] = corr(IMUpost2,Qpost2);
OUT.Corr_speed_FR_r_p_post2 = [rpost2; ppost2];
[rpost3, ppost3] = corr(IMUpost3,Qpost3);
OUT.Corr_speed_FR_r_p_post3 = [rpost3; ppost3];


if PLOT_IT
    
    figure
    subplot(4,1,1)
    plot(IMUall_x_uS(:,1)/60e6,IMUall)
    plot_markers_simple(TIMES.EventStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    xlabel('min')
    ylabel('IMU speed')
    title(sprintf('Binsize %d s',binsize_Q_S_ms/1000))
    
    subplot(4,1,2)
    imagesc(start_end_times(:,1)/60e6,[],Q')
    plot_markers(TIMES.EventStartEndUsec/60e6);
    colorbar
    xlabel('min')
    ylabel('neuron')
    title(sprintf('Binsize %d s',binsize_Q_S_ms/1000))
    plot_markers_simple(TIMES.EventStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    subplot(4,1,3)
    % same but subtract mean rate at beginning.
    
    imagesc(start_end_times(:,1)/60e6,[],Qn')
    plot_markers_simple(TIMES.EventStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    colorbar
    xlabel('min')
    ylabel('neuron')
    title('Subtract rate before injection.')
    
    subplot(4,1,4)
    plot_confidence_intervals(start_end_times(:,1)/60e6,Q')
    ylabel('Mean Neural activity')
    yyaxis right
    plot_confidence_intervals(start_end_times(:,1)/60e6,Qn',[],[.8 .1 .1])
    plot_markers_simple(TIMES.EventStartEndUsec/60e6);
    ylabel('norm neural activity')
end

% plot the firing and movement
if PLOT_IT
    clrs = lines(5);
    figure
    plot_confidence_intervals(Qall_x_uS(:,1)/60e6,Qall',[],clrs(1,:))
    yyaxis right
    plot(IMUall)
    plot_markers_simple(TIMES.EventStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    % plot just the ketamine period
    figure
    plot_confidence_intervals(x2(:,1)/60e6,Qpost2',[],clrs(2,:))
    yyaxis right
    plot(y2(:,1)/60e6,IMUpost2)

end

%% Get autocorrelation during different periods for neuron classification

[AC,x] = AutoCorrArray_Abhi(TS,bin_size_ac_ms*1000,acorr_width_ms*1000,[TIMES.BaselineUsec;TIMES.PostInjectionUsec]);
OUT.AC_base_post1_post2_post3 = AC;
OUT.AC_x_ms = x/1000;
AC_x_ms = x/1000;
% plot the AC
% 
% figure
% imagesc(squeeze(norm_AC(:,:,1)))

for iv = 1:4
    ix = find(AC_x_ms > max(AC_x_ms)-100);
    norm_AC(:,:,iv) = AC(:,:,iv) - repmat(mean(AC(:,ix,iv),2),1,Cols(AC(:,:,iv))); titstr = 'Norm By Mean Rate';
    norm_AC(:,:,iv) = AC(:,:,iv) ./ repmat(sum(AC(:,:,iv),2),1,Cols(AC(:,:,iv)));
    %     norm_AC(:,:,iv) = standardize_range(norm_AC(:,:,iv)')';
end

OUT.AC_norm_base_post123 = norm_AC;

%% Look at local variance...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('LV')
TSx{1} = Restrict(TS,TIMES.BaselineUsec);
TSx{2} = Restrict(TS,TIMES.PostInjectionUsec(1,:));
TSx{3} = Restrict(TS,TIMES.PostInjectionUsec(2,:));
TSx{4} = Restrict(TS,TIMES.PostInjectionUsec(3,:));

SSS = [];
OUT.LocVar_base_post1_post2_post3 = [];
for ii = 1:length(TSx{1})
    for jj = 1:length(TSx)
        OUT.LocVar_base_post1_post2_post3(ii,jj) = LocalVariance(diff(TSx{jj}{ii})/1000,1e4); % ISI is in ms and outlier threshold is 10 sec
        SSS{ii,jj} = Spiking_summary_stats(TSx{jj}{ii}/1000);
    end
end
OUT.SSS = SSS;

% Smooth the AC.
% norm_smth_AC = sgolayfilt(norm_AC',3,9)';

% [mn_AC, ci_AC]= normci(squeeze(norm_AC(:,:,1)));
% [mn_AC3, ci_AC3]= normci(squeeze(norm_AC(:,:,3)));
% 
% figure
% plot_confidence_intervals(AC_x_ms, mn_AC, ci_AC, clrs(1,:))
% hold on
% plot_confidence_intervals(AC_x_ms, mn_AC3, ci_AC3, clrs(3,:))
% xlabel('time (ms)')
% ylabel('AC')
% title('Rat 320 Ses 24 Auto Corr for all neurons during baseline and ketamine')
% pubify_figure_axis

%% testing how to plot the movement and spike activity 

% [rho, pval] = corr(IMUall,Qall);
% figure; scatter(IMUall,Qall(:,5)) % plotting the ketamine period IMU and one neuron activity
% hold on
% P = polyfit(IMUpost2, Qpost2(:,5), 1); % getting the linear fit
% Bfit = polyval(P, IMUpost2); % creating the line?? pulled this from the interwebs
% plot(IMUpost2, Bfit, '-r')
% 
% 
% figure; scatter(IMUpost2,Qpost2(:,5)) % plotting the ketamine period IMU and one neuron activity
% hold on
% P = polyfit(IMUpost2, Qpost2(:,5), 1); % getting the linear fit
% Bfit = polyval(P, IMUpost2); % creating the line?? pulled this from the interwebs
% plot(IMUpost2, Bfit, '-r')
% ylabel('Firirng rate .5s bins')
% xlabel('IMU speed .5s bins')
% title('Rat 342 Session 4 Neuron 5 ketamine period movement and spike activity, r = .24 p = <.001')
% pubify_figure_axis
% 
% [rho1, pval1] = corr(IMUpost1,Qpost1);
% 
% figure; scatter(IMUpost1,Qpost1(:,30)) % plotting the ketamine period IMU and one neuron activity
% hold on
% P = polyfit(IMUpost1,Qpost1(:,30), 1); % getting the linear fit
% Bfit = polyval(P, IMUpost1); % creating the line?? pulled this from the interwebs
% plot(IMUpost1, Bfit, '-r')
% ylabel('Firirng rate .5s bins')
% xlabel('IMU speed .5s bins')
% title('Rat 342 Session 4 Neuron 30 peak dyskinetic period movement and spike activity, r = .23 p = <.001')
% pubify_figure_axis
% 
% % PCA analysis for neurons in one session
% [coeff,score,latent,tsquared,explained] = pca(Qpost2);
% 
% figure
% plot(explained)
% xlabel('PC')
% ylabel('Percent Total variance Explained')
% title('Rat 342 Dyskinetic period')
% 
% figure
% subplot(3,1,1)
% plot(score(:,1))
% title(sprintf('PC1 %d variance explained',explained(1)))
% subplot(3,1,2)
% plot(score(:,2))
% title(sprintf('PC2 %d variance explained',explained(2)))
% subplot(3,1,3)
% plot(score(:,3))
% title(sprintf('PC3 %d variance explained',explained(3)))
% xlabel('PostKetamine time')
% 
% [rho,pval] = corr(IMUpost1,score(:,2));
% 
% figure; scatter(IMUpost1,score(:,1)) 
% hold on
% xlabel('IMUspeed')
% ylabel('PC1 of all Ses 4 neurons')
% P = polyfit(IMUpost1, score(:,1), 1);
% Bfit = polyval(P, IMUpost1);
% plot(IMUpost1, Bfit, '-r')
% title('PC1 Rat 342 Ses 4 neurons Dyskinetic period, r = -.09 p <.001')
% pubify_figure_axis
% 
% figure; scatter(IMUpost1,score(:,2)) 
% hold on
% xlabel('IMUspeed')
% ylabel('PC2 of all Ses 4 neurons')
% P = polyfit(IMUpost1, score(:,2), 1);
% Bfit = polyval(P, IMUpost1);
% plot(IMUpost1, Bfit, '-r')
% title('PC2 Rat 342 Ses 4 neurons Dyskinetic period, r = .06 p <.001')
% pubify_figure_axis
% 
% figure; scatter(IMUpost2,score(:,3)) 
% hold on
% xlabel('IMUspeed')
% ylabel('PC3 of all Ses 3 neurons')
% P = polyfit(IMUpost2, score(:,3), 1);
% Bfit = polyval(P, IMUpost2);
% plot(IMUpost2, Bfit, '-r')
% title('PC3 Rat 350 Ses 3 neurons, r = -.069 p <.001')
% 
% % plotting eigen vectors
% figure
% bar(coeff(:,1))
% 
% % plotting the PC for the entire recording session
% 
% [coeff,score,latent,tsquared,explained] = pca(Qall);
% 
% figure
% scatter(score(:,1),score(:,2))
% 
% % Creating groups
% test = Qall_x_uS(:,1) > TIMES.BaselineUsec(1) & Qall_x_uS(:,1) <= TIMES.BaselineUsec(2);
% test1 = Qall_x_uS(:,1) > TIMES.PostInjectionUsec(1,1) & Qall_x_uS(:,1) <= TIMES.PostInjectionUsec(1,2);
% test2 = Qall_x_uS(:,1) > TIMES.PostInjectionUsec(2,1) & Qall_x_uS(:,1) <= TIMES.PostInjectionUsec(2,2);
% test3 = Qall_x_uS(:,1) > TIMES.PostInjectionUsec(3,1) & Qall_x_uS(:,1) <= TIMES.PostInjectionUsec(3,2);
% 
% test = double(test);
% test(test1) = 2;
% test(test2) = 3;
% test(test3) = 4;
% drug_g = test;
% 
% clrs = lines(5);
% figure
% gscatter(score(:,1),score(:,2),drug_g,clrs(1:5,:),'.o+x*',5)
% xlabel('PC1')
% ylabel('PC2')
% title('Rat 320 Ses 24 Entire session PCA')
% pubify_figure_axis
% % Correlation before and after ketamine 
% R_base = corrcoef(Qbase);
% figure 
% imagesc(R_base)
% colorbar
% caxis(prctile(R_base(:),[1 99]))
% 
% R_ldo = corrcoef(Qpost1);
% figure 
% imagesc(R_ldo)
% colorbar
% caxis(prctile(R_ldo(:),[1 99]))
% 
% R_ket = corrcoef(Qpost2);
% figure 
% imagesc(R_ket)
% colorbar
% caxis(prctile(R_ket(:),[1 99]))
