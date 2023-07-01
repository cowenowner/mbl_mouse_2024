function OUT = Q2_ketamine_m1_vs_striatum_single_unit()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is there a significant difference in firing rate of neurons in Motor
% cortex and Striatum during pre and post ketamine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = [];
bin_size_ms = 2;
acorr_width_ms = 200;       % defining bin sizes and things
PLOT_IT = true;
intervals_around_ket_min = [-20 -5; 5 20; 25 40];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function loads a ton of things that are typically used
% for any analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try
    [GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;
% catch
%     disp('Failed to load files on first try for some reason. Box?')
%     pause(10)
%     [GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;
% end
    
OUT.DEPTHS = DEPTHS; % Send this info out of the function for meta analysis.
OUT.META = META;
OUT.EVT = EVT;
OUT.SP = SP;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Seperating the neurons in Striatum and Motor cortex
SPT = struct2table(SP);
SIX = SPT.Depth_uM > 3000;
MIX = SPT.Depth_uM < 3000;
TS_STR = TS(:,SIX);
TS_M1 = TS(:,MIX);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the relevant times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TIMES.KetStartEndUsec(1,1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
TIMES.KetStartEndUsec(1,2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
TIMES.BaselineUsec(1,1) = TIMES.KetStartEndUsec(1) + intervals_around_ket_min(1,1)*60*1e6;
TIMES.BaselineUsec(1,2) = TIMES.KetStartEndUsec(1) + intervals_around_ket_min(1,2)*60*1e6;
TIMES.PostInjectionUsec(1,1) = TIMES.KetStartEndUsec(2) + intervals_around_ket_min(2,1)*60*1e6;
TIMES.PostInjectionUsec(1,2) = TIMES.KetStartEndUsec(2) + intervals_around_ket_min(2,2)*60*1e6;
TIMES.PostInjectionUsec(2,1) = TIMES.KetStartEndUsec(2) + intervals_around_ket_min(3,1)*60*1e6;
TIMES.PostInjectionUsec(2,2) = TIMES.KetStartEndUsec(2) + intervals_around_ket_min(3,2)*60*1e6;

OUT.TM = TIMES;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load inertial, POS, and calculate string pull rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMU = LK_Load_and_Process_IMU('Inertial_data.mat');
% POS = LK_Load_and_Clean_POS('POS.mat');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Find indices for the times for baseline, and the two post-ketamine periods.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IX1 = POS.Time_uS > TIMES.BaselineUsec(1) & POS.Time_uS < TIMES.BaselineUsec(2);
% IX2 = POS.Time_uS > TIMES.PostInjectionUsec(1,1) & POS.Time_uS < TIMES.PostInjectionUsec(1,2);
% IX3 = POS.Time_uS > TIMES.PostInjectionUsec(2,1) & POS.Time_uS < TIMES.PostInjectionUsec(2,2);
% OUT.Speed_Before_Early_Late = [nanmean(POS.speed(IX1)) nanmean(POS.speed(IX2)) nanmean(POS.speed(IX3))];
% 
% IX1 = IMU.t_uS > TIMES.BaselineUsec(1) & IMU.t_uS < TIMES.BaselineUsec(2);
% IX2 = IMU.t_uS > TIMES.PostInjectionUsec(1,1) & IMU.t_uS < TIMES.PostInjectionUsec(1,2);
% IX3 = IMU.t_uS > TIMES.PostInjectionUsec(2,1) & IMU.t_uS < TIMES.PostInjectionUsec(2,2);
% % Store information for analysis
% OUT.Jerk_Before_Early_Late = [nanmean(IMU.absjerk(IX1)) nanmean(IMU.absjerk(IX2)) nanmean(IMU.absjerk(IX3))];
% OUT.JerkPC1_Before_Early_Late = [nanmean(IMU.absjerkpc(IX1)) nanmean(IMU.absjerkpc(IX2)) nanmean(IMU.absjerkpc(IX3))];
% OUT.IMU_Speed_Before_Early_Late = [nanmean(IMU.speed(IX1)) nanmean(IMU.speed(IX2)) nanmean(IMU.speed(IX3))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Always visualize.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movement
% if PLOT_IT
%     figure
%     plot(IMU.t_uS/60e6,IMU.absjerk)
%     hold on
%     plot(IMU.t_uS/60e6,IMU.speed,'k')
%     
%     ylabel('Jerk')
%     xlabel('min')
%     plot_markers_simple(TIMES.KetStartEndUsec/60e6,[],[],GP.Colors.Ketamine)
%     plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'g')
%     plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'c')
%     
%     yyaxis right
%     plot(POS.Time_uS/60e6,POS.speed)
%     ylabel('speed from x and y')
%     pubify_figure_axis
%     axis tight
%     title('Position Data')
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % speed
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure
%     subplot(3,1,1)
%     bar(OUT.Speed_Before_Early_Late)
%     set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
%     ylabel('Speed')
%     pubify_figure_axis
%     subplot(3,1,2)
%     bar(OUT.Jerk_Before_Early_Late)
%     set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
%     ylabel('Jerk')
%     pubify_figure_axis
%     subplot(3,1,3)
%     bar(OUT.IMU_Speed_Before_Early_Late)
%     set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
%     ylabel('IMU Speed')
%     pubify_figure_axis
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overview: Generate a neuron by time matrix and mark when ketamine was delivered.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binsize_ms = 10000;
bins_s = binsize_ms/1000;

% Striatum
% [Q_S,start_end_times_S] = Bin_ts_array(TS(:,SIX), binsize_ms*1000); % workhorse function for neuron x time matrix
% IX = start_end_times_S(:,1)<TIMES.KetStartEndUsec(1);
% Qn_S = Q_S-mean(Q_S(IX,:));
% [~,sc_S] = pca(Qn_S);
% 
% edges = TIMES.BaselineUsec(1,1):binsize_ms*1000:TIMES.BaselineUsec(1,2);
% [Qbase_S,bins_S_uS] = Bin_ts_array(TS(:,SIX), edges);
% edges = TIMES.PostInjectionUsec(1,1):binsize_ms*1000:TIMES.PostInjectionUsec(1,2);
% [Qpost1_S,x1_S] = Bin_ts_array(TS(:,SIX), edges);
% edges = TIMES.PostInjectionUsec(2,1):binsize_ms*1000:TIMES.PostInjectionUsec(2,2);
% [Qpost2_S,x2_S] = Bin_ts_array(TS(:,SIX), edges);

[Q,start_end_times] = Bin_ts_array(TS_STR, binsize_ms*1000); % workhorse function for neuron x time matrix
IX = start_end_times(:,1)<TIMES.KetStartEndUsec(1);
Qn = Q-mean(Q(IX,:));
[~,sc] = pca(Qn);

edges = TIMES.BaselineUsec(1,1):binsize_ms*1000:TIMES.BaselineUsec(1,2);
[Qbase,bins_uS] = Bin_ts_array(TS_STR, edges);
edges = TIMES.PostInjectionUsec(1,1):binsize_ms*1000:TIMES.PostInjectionUsec(1,2);
[Qpost1,x1] = Bin_ts_array(TS_STR, edges);
edges = TIMES.PostInjectionUsec(2,1):binsize_ms*1000:TIMES.PostInjectionUsec(2,2);
[Qpost2,x2] = Bin_ts_array(TS_STR, edges);

mbase = mean(Qbase/bins_s);
m1 = mean(Qpost1/bins_s);
m2 = mean(Qpost2/bins_s);

OUT.Striatum.FRates_base = mbase;
OUT.Striatum.FRates_post1 = m1;
OUT.Striatum.FRates_post2 = m2;

OUT.Striatum.Qbase = Qbase;
OUT.Striatum.Qbase_bins_uS = bins_uS;
OUT.Striatum.Qpost1 = Qpost1;
OUT.Striatum.Qpost1_bins_uS = x1;
OUT.Striatum.Qpost2 = Qpost2;
OUT.Striatum.Qpost2_bins_uS = x2;

if PLOT_IT

    figure
    subplot(3,1,1)
    imagesc(start_end_times(:,1)/60e6,[],Q')
    plot_markers(TIMES.KetStartEndUsec/60e6);
    colorbar
    xlabel('min')
    ylabel('neuron')
    title(sprintf('Striatum; Binsize %d s',binsize_ms/1000))
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    subplot(3,1,2)
    % same but subtract mean rate at beginning.
    
    imagesc(start_end_times(:,1)/60e6,[],Qn')
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    colorbar
    xlabel('min')
    ylabel('neuron')
    title('Striatum; Subtract rate before injection.')
    
    subplot(3,1,3)
    plot_confidence_intervals(start_end_times(:,1)/60e6,Q')
    ylabel('Mean Neural activity')
    yyaxis right
    plot_confidence_intervals(start_end_times(:,1)/60e6,Qn',[],[.8 .1 .1])
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    ylabel('norm neural activity')
    title('Striatum neural activity')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot mean firing rates before and after ketamine.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    subplot(2,1,1)
    
    %     plot([mbase;m1;m2]')
    plot(start_end_times(:,1)/60e6,sc(:,1:2));
    axis tight
    %     ylabel('Hz')
    %     legend('pre','p1','p2','FontSize',8)
    %     title('Firing rates')
    subplot(2,1,2)
    h = bar(([m1;m2]- mbase)');
    h(1).FaceColor = 'k';
    h(2).FaceColor = 'b';
    set(gca,'Colormap',jet())
    legend('p1','p2','FontSize',8)
    legend boxoff
    title('Striatum; rate after subtracting baseline')
    xlabel('Neuron')
    
end

% Motor cortex
% [Q_M,start_end_times_M] = Bin_ts_array(TS(:,MIX), binsize_ms*1000); % workhorse function for neuron x time matrix
% IX = start_end_times_M(:,1)<TIMES.KetStartEndUsec(1);
% Qn_M = Q_M-mean(Q_M(IX,:));
% [~,sc_M] = pca(Qn_M);
% 
% edges = TIMES.BaselineUsec(1,1):binsize_ms*1000:TIMES.BaselineUsec(1,2);
% [Qbase_M,bins_M_uS] = Bin_ts_array(TS(:,MIX), edges);
% edges = TIMES.PostInjectionUsec(1,1):binsize_ms*1000:TIMES.PostInjectionUsec(1,2);
% [Qpost1_M,x1_M] = Bin_ts_array(TS(:,MIX), edges);
% edges = TIMES.PostInjectionUsec(2,1):binsize_ms*1000:TIMES.PostInjectionUsec(2,2);
% [Qpost2_M,x2_M] = Bin_ts_array(TS(:,MIX), edges);
% 
% mbase_M = mean(Qbase_M/bins_s);
% m1_M = mean(Qpost1_M/bins_s);
% m2_M = mean(Qpost2_M/bins_s);

[Q,start_end_times] = Bin_ts_array(TS_M1, binsize_ms*1000); % workhorse function for neuron x time matrix
IX = start_end_times(:,1)<TIMES.KetStartEndUsec(1);
Qn = Q-mean(Q(IX,:));
[~,sc] = pca(Qn);

edges = TIMES.BaselineUsec(1,1):binsize_ms*1000:TIMES.BaselineUsec(1,2);
[Qbase,bins_uS] = Bin_ts_array(TS_M1, edges);
edges = TIMES.PostInjectionUsec(1,1):binsize_ms*1000:TIMES.PostInjectionUsec(1,2);
[Qpost1,x1] = Bin_ts_array(TS_M1, edges);
edges = TIMES.PostInjectionUsec(2,1):binsize_ms*1000:TIMES.PostInjectionUsec(2,2);
[Qpost2,x2] = Bin_ts_array(TS_M1, edges);

mbase = mean(Qbase/bins_s);
m1 = mean(Qpost1/bins_s);
m2 = mean(Qpost2/bins_s);

OUT.Motor.FRates_base = mbase;
OUT.Motor.FRates_post1 = m1;
OUT.Motor.FRates_post2 = m2;

OUT.Motor.Qbase = Qbase;
OUT.Motor.Qbase_bins_uS = bins_uS;
OUT.Motor.Qpost1 = Qpost1;
OUT.Motor.Qpost1_bins_uS = x1;
OUT.Motor.Qpost2 = Qpost2;
OUT.Motor.Qpost2_bins_uS = x2;

if PLOT_IT
    
    figure
    subplot(3,1,1)
    imagesc(start_end_times(:,1)/60e6,[],Q')
    plot_markers(TIMES.KetStartEndUsec/60e6);
    colorbar
    xlabel('min')
    ylabel('neuron')
    title(sprintf('Motor Cortex; Binsize %d s',binsize_ms/1000))
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    
    subplot(3,1,2)
    % same but subtract mean rate at beginning.
    
    imagesc(start_end_times(:,1)/60e6,[],Qn')
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    colorbar
    xlabel('min')
    ylabel('neuron')
    title('Motor cortex; Subtract rate before injection.')
    
    subplot(3,1,3)
    plot_confidence_intervals(start_end_times(:,1)/60e6,Q')
    ylabel('Mean Neural activity')
    yyaxis right
    plot_confidence_intervals(start_end_times(:,1)/60e6,Qn',[],[.8 .1 .1])
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    ylabel('norm neural activity')
    title('Motor cortex neural activity')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot mean firing rates before and after ketamine.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure
    subplot(2,1,1)
    
    %     plot([mbase;m1;m2]')
    plot(start_end_times(:,1)/60e6,sc(:,1:2));
    axis tight
    %     ylabel('Hz')
    %     legend('pre','p1','p2','FontSize',8)
    %     title('Firing rates')
    subplot(2,1,2)
    h = bar(([m1;m2]- mbase)');
    h(1).FaceColor = 'k';
    h(2).FaceColor = 'b';
    set(gca,'Colormap',jet())
    legend('p1','p2','FontSize',8)
    legend boxoff
    title('Motor cortex; rate after subtracting baseline')
    xlabel('Neuron')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot overall difference score in firing rates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure
%     histogram_cowen({[m1 - mbase], [m2 - mbase]},.1);
%     legend('p1','p2','FontSize',8)
%     legend boxoff
%     plot_vert_line_at_zero
%     xlabel('Change in Firing Rate')
%     [~,p1] = ttest([m1 - mbase]);
%     [~,p2] = ttest([m2 - mbase]);
%     title(sprintf('ttest: p1 %0.4f, p2 %0.4f',p1,p2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at autocorrelations before and after and the difference between
% do it for the pre-injection period but also for the 2-20 minute and 30-50
% minute post periods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[AC,x] = AutoCorrArray(TS_STR,bin_size_ms*1000,acorr_width_ms*1000,[TIMES.BaselineUsec;TIMES.PostInjectionUsec]);
ACd = AC -  AC(:,:,1);
cm = lines(5);
OUT.Striatum.AC = AC;
OUT.Striatum.AC_x_ms = x/1000;
OUT.AC_conditions = {'pre' 'post 1' 'post 2'};
% disp('ACorr')
if PLOT_IT
    titls = OUT.AC_conditions;
    CIX = x/1000 > 120;
    
    figure
    for ii = 1:length(titls)
        subplot(1,length(titls),ii)
        a = squeeze(AC(:,:,ii));
        a = a - mean(a(:,CIX),2);
        imagesc(x/1000,[],a)
        title(titls{ii})
    end
    equalize_color_axes
    
    figure
    titls = {'post 1' 'post 2'};
    for ii = 1:length(titls)
        hh(ii) = subplot(1,length(titls),ii);
        imagesc(x/1000,[],AC(:,:,ii+1)- AC(:,:,1))
        title(titls{ii})
        colorbar
        colormap(parula);
    end
    equalize_color_axes(hh);
    
    figure
    for ii = 1:3
        a = squeeze(AC(:,:,ii));
        a = a - mean(a(:,CIX),2);
        plot_confidence_intervals(x/1000,a,[],cm(ii,:));
    end
    legend_color_text(OUT.AC_conditions,cm(1:3,:));
    
    figure
    plot_confidence_intervals(x/1000,squeeze(AC(:,:,1)),[],cm(1,:))
    plot_confidence_intervals(x/1000,squeeze(AC(:,:,3)),[],cm(2,:))
    plot_horiz_line_at_zero
    
    
    figure
    plot_confidence_intervals(x/1000,squeeze(ACd(:,:,2)),[],cm(1,:))
    plot_confidence_intervals(x/1000,squeeze(ACd(:,:,3)),[],cm(2,:))
    plot_horiz_line_at_zero
    xlabel('ms')
    ylabel('Difference from baseline')
    title('Straitum; Autocorr after ketamine relative to baseline')
    legend_color_text({'p1' 'p2'},{cm(1,:) cm(2,:)});
end

%%

[AC,x] = AutoCorrArray(TS_M1,bin_size_ms*1000,acorr_width_ms*1000,[TIMES.BaselineUsec;TIMES.PostInjectionUsec]);
ACd = AC -  AC(:,:,1);
cm = lines(5);
OUT.Motor.AC = AC;
OUT.Motor.AC_x_ms = x/1000;

% disp('ACorr')
if PLOT_IT
    titls = OUT.AC_conditions;
    CIX = x/1000 > 120;
    
    figure
    for ii = 1:length(titls)
        subplot(1,length(titls),ii)
        a = squeeze(AC(:,:,ii));
        a = a - mean(a(:,CIX),2);
        imagesc(x/1000,[],a)
        title(titls{ii})
    end
    equalize_color_axes
    
    figure
    titls = {'post 1' 'post 2'};
    for ii = 1:length(titls)
        hh(ii) = subplot(1,length(titls),ii);
        imagesc(x/1000,[],AC(:,:,ii+1)- AC(:,:,1))
        title(titls{ii})
        colorbar
        colormap(parula);
    end
    equalize_color_axes(hh);
    
    figure
    for ii = 1:3
        a = squeeze(AC(:,:,ii));
        a = a - mean(a(:,CIX),2);
        plot_confidence_intervals(x/1000,a,[],cm(ii,:));
    end
%     legend_color_text(OUT.AC_conditions,cm(1:3,:));
    
    figure
    plot_confidence_intervals(x/1000,squeeze(AC(:,:,1)),[],cm(1,:))
    plot_confidence_intervals(x/1000,squeeze(AC(:,:,3)),[],cm(2,:))
    plot_horiz_line_at_zero
    
    
    figure
    plot_confidence_intervals(x/1000,squeeze(ACd(:,:,2)),[],cm(1,:))
    plot_confidence_intervals(x/1000,squeeze(ACd(:,:,3)),[],cm(2,:))
    plot_horiz_line_at_zero
    xlabel('ms')
    ylabel('Difference from baseline')
    title('Motor cortex; Autocorr after ketamine relative to baseline')
    legend_color_text({'p1' 'p2'},{cm(1,:) cm(2,:)});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smooth autocorr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('SACorr')

% binsize_ms = 2;
% smth_win_ms = 8;
% acorr_size_ms = 250;
% jitter_ms = 500;
% if 0
%     for ii = 1:length(TS)
%         for jj = 1:length(TSx)
%             [OUT.AcorrSmth_Before_Early_Late(ii,:,jj),OUT.AcorrSmth_x_ms, CI] = ...
%                 CrossCorr_smooth(TSx{jj}{ii}/1000,[],binsize_ms,acorr_size_ms/binsize_ms,smth_win_ms/binsize_ms,jitter_ms);
%             OUT.AcorrSmthCI95up_Before_Early_Late(ii,:,jj)= CI(1,:);
%             OUT.AcorrSmthCI95down_Before_Early_Late(ii,:,jj)= CI(2,:);
%             %         CrossCorr_smooth(TSx{jj}{ii}/1000,[],binsize_ms,acorr_size_ms/binsize_ms,smth_win_ms/binsize_ms,jitter_ms)
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Look at local variance...
% %%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('LV')

TSx_STR{1} = Restrict(TS_STR,TIMES.BaselineUsec);
TSx_STR{2} = Restrict(TS_STR,TIMES.PostInjectionUsec(1,:));
TSx_STR{3} = Restrict(TS_STR,TIMES.PostInjectionUsec(2,:));

SSS = [];
OUT.Striatum.LocVar_Before_Early_Late = [];
for ii = 1:length(TSx_STR{1})
    for jj = 1:length(TSx_STR)
        OUT.Striatum.LocVar_Before_Early_Late(ii,jj) = LocalVariance(diff(TSx_STR{jj}{ii}),10e6);
        SSS{ii,jj} = Spiking_summary_stats(TSx_STR{jj}{ii}/1000);
    end
end
OUT.Striatum.SSS = SSS;

TSx_M1{1} = Restrict(TS_M1,TIMES.BaselineUsec);
TSx_M1{2} = Restrict(TS_M1,TIMES.PostInjectionUsec(1,:));
TSx_M1{3} = Restrict(TS_M1,TIMES.PostInjectionUsec(2,:));

SSS = [];
OUT.Motor.LocVar_Before_Early_Late = [];
for ii = 1:length(TSx_M1{1})
    for jj = 1:length(TSx_M1)
        OUT.Motor.LocVar_Before_Early_Late(ii,jj) = LocalVariance(diff(TSx_M1{jj}{ii}),10e6);
        SSS{ii,jj} = Spiking_summary_stats(TSx_M1{jj}{ii}/1000);
    end
end
OUT.Motor.SSS = SSS;
