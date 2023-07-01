function OUT = Q1_What_does_ketamine_do_firing_rate_Abhi()
% function OUT = Q1_What_does_ketamine_do_firing_rate()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How does ketamine change the firing rate properties and other
% spike-statistics of neurons? Autocorrs, etc...
%
% Other functions will deal with higher-order statistics like ensemble
% activity or LFP coupling. We will keep this one restricted to keep things
% simple.
%
% Determine if single-units oscillate (autocorr) more after ketamine injection.
% Compare autocorrs before and after injection.
%
% NOTE: All Q functions must be able to run in each Session directory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if nargin < 1
%     type = 'ketamine';
% end
OUT = [];
bin_size_ac_ms = 2;
% binsize_small_ms = 100;
% binsize_Q_ms = 10000; % 10 sec for FR change analysis
binsize_Q_ms = 20; % 20 msec for FR change analysis
% SPIKE_PSD = false; % good but takes A LONG time. run this only overnight.
acorr_width_ms = 200;       % defining bin sizes and things
PLOT_IT = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function loads a ton of things that are typically used
% for any analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;
OUT.SP = SP;
% catch
%     disp('Failed to load files on first try for some reason. Box?')
%     pause(10)
%     [GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;
% end
OUT.binsize_Q_ms = binsize_Q_ms;
OUT.bin_size_ac_ms = bin_size_ac_ms;

OUT.DEPTHS = DEPTHS; % Send this info out of the function for meta analysis.
OUT.META = META;
OUT.EVT = EVT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the relevant times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inj1_start_end_uS = [];
inj2_start_end_uS = [];

if any(E.EventID == 'KetInjectionStart') && any(E.EventID == 'SalineInjectionStart')
        inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
        inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
        inj1_start_end_uS(1) = E.MinFromStart(E.EventID == 'SalineInjectionStart')*60*1e6;
        inj1_start_end_uS(2) = E.MinFromStart(E.EventID == 'SalineInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 2 22; 60 80];
        big_peri_event_min = [-25 110];
        OUT.inj1 = 'Sal';
        OUT.inj2 = 'Ket';

elseif any(E.EventID == 'LDOPAInjectionStart')
        inj1_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
        inj1_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        if any(E.EventID == 'SalineInjectionStart')
            inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'SalineInjectionStart')*60*1e6;
            inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'SalineInjectionEnd')*60*1e6;
            OUT.inj2 = 'Sal';
        elseif any(E.EventID == 'KetInjectionStart')
            inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
            inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
            OUT.inj2 = 'Ket';
        end
        intervals_around_evt_min =  [-25 -5; 2 22; 60 80; 5 25];
        big_peri_event_min = [-25 180];
        OUT.inj1 = 'Ldo';
end

% return
OUT.intervals_around_evt_min = intervals_around_evt_min;
TIMES.Event1StartEndUsec(1,1) = inj1_start_end_uS(1);
TIMES.Event1StartEndUsec(1,2) = inj1_start_end_uS(2);
TIMES.Event2StartEndUsec(1,1) = inj2_start_end_uS(1);
TIMES.Event2StartEndUsec(1,2) = inj2_start_end_uS(2);
TIMES.BaselineUsec(1,1) = TIMES.Event1StartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.BaselineUsec(1,2) = TIMES.Event1StartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.PostInjectionUsec(1,1) = TIMES.Event2StartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.PostInjectionUsec(1,2) = TIMES.Event2StartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.PostInjectionUsec(2,1) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(2,1)*60*1e6;
TIMES.PostInjectionUsec(2,2) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(2,2)*60*1e6;
TIMES.PostInjectionUsec(3,1) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(3,1)*60*1e6;
TIMES.PostInjectionUsec(3,2) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(3,2)*60*1e6;


TIMES.PeriKetamineUsec = [TIMES.Event1StartEndUsec(1) + big_peri_event_min(1)*60*1e6 TIMES.Event1StartEndUsec(1) + big_peri_event_min(2)*60*1e6 ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load inertial, POS, and calculate string pull rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IMU = LK_Load_and_Process_IMU('Inertial_data.mat');
POS = LK_Load_and_Clean_POS_Abhi('POS.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find indices for the times for baseline, and the two post-ketamine periods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ix1 = POS.Time_uS > TIMES.BaselineUsec(1) & POS.Time_uS < TIMES.BaselineUsec(2);
ix2 = POS.Time_uS > TIMES.PostInjectionUsec(1,1) & POS.Time_uS < TIMES.PostInjectionUsec(1,2);
ix3 = POS.Time_uS > TIMES.PostInjectionUsec(2,1) & POS.Time_uS < TIMES.PostInjectionUsec(2,2);
ix4 = POS.Time_uS > TIMES.PostInjectionUsec(3,1) & POS.Time_uS < TIMES.PostInjectionUsec(3,2);
OUT.Speed_Base_Pre_Early_Late = [nanmean(POS.speed(ix1)) nanmean(POS.speed(ix2)) nanmean(POS.speed(ix3)) nanmean(POS.speed(ix4))];

IX1 = IMU.t_uS > TIMES.BaselineUsec(1) & IMU.t_uS < TIMES.BaselineUsec(2);
IX2 = IMU.t_uS > TIMES.PostInjectionUsec(1,1) & IMU.t_uS < TIMES.PostInjectionUsec(1,2);
IX3 = IMU.t_uS > TIMES.PostInjectionUsec(2,1) & IMU.t_uS < TIMES.PostInjectionUsec(2,2);
IX4 = IMU.t_uS > TIMES.PostInjectionUsec(3,1) & IMU.t_uS < TIMES.PostInjectionUsec(3,2);

% Store information for analysis
OUT.Jerk_base_pre_early_late = [nanmean(IMU.absjerk(IX1)) nanmean(IMU.absjerk(IX2)) nanmean(IMU.absjerk(IX3)) nanmean(IMU.absjerk(IX4))];
OUT.JerkPC1_base_pre_early_late = [nanmean(IMU.absjerkpc(IX1)) nanmean(IMU.absjerkpc(IX2)) nanmean(IMU.absjerkpc(IX3)) nanmean(IMU.absjerkpc(IX4))];
OUT.IMU_Speed_base_pre_early_late = [nanmean(IMU.speed(IX1)) nanmean(IMU.speed(IX2)) nanmean(IMU.speed(IX3)) nanmean(IMU.speed(IX4))];

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
    subplot(3,1,1)
    bar(OUT.Speed_Base_Sal_Early_Late)
    set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
    ylabel('Speed')
    pubify_figure_axis
    subplot(3,1,2)
    bar(OUT.Jerk_base_post1_post2_post3)
    set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
    ylabel('Jerk')
    pubify_figure_axis
    subplot(3,1,3)
    bar(OUT.IMU_Speed_base_post1_post2_post3)
    set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
    ylabel('IMU Speed')
    pubify_figure_axis
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overview: Generate a neuron by time matrix and mark when ketamine was delivered.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Q,start_end_times] = Bin_ts_array(TS, binsize_Q_ms*1000); % workhorse function for neuron x time matrix
IX = start_end_times(:,1)<TIMES.Event1StartEndUsec(1);
Qn = Q-mean(Q(IX,:));
[~,sc] = pca(Qn);

edges = TIMES.BaselineUsec(1,1):binsize_Q_ms*1000:TIMES.BaselineUsec(1,2);
[Qbase,bins_uS] = Bin_ts_array(TS, edges);
edges = TIMES.PostInjectionUsec(1,1):binsize_Q_ms*1000:TIMES.PostInjectionUsec(1,2);
[Qpost1,x1] = Bin_ts_array(TS, edges);
edges = TIMES.PostInjectionUsec(2,1):binsize_Q_ms*1000:TIMES.PostInjectionUsec(2,2);
[Qpost2,x2] = Bin_ts_array(TS, edges);
edges = TIMES.PostInjectionUsec(3,1):binsize_Q_ms*1000:TIMES.PostInjectionUsec(3,2);
[Qpost3,x3] = Bin_ts_array(TS, edges);

edges = TIMES.PeriKetamineUsec(1):binsize_Q_ms*1000:TIMES.PeriKetamineUsec(2);
[Qall,Qall_x_uS] = Bin_ts_array(TS, edges);
Qall = int16(Qall);
TSaligned = cell(size(TS));
for ii = 1:length(TS)
    if ~isempty(TS{ii})
        TSaligned{ii} = TS{ii} - TIMES.Event2StartEndUsec(2);
    end
end
OUT.TSaligned = TSaligned;

bins_s = binsize_Q_ms/1000;
mbase = mean(Qbase/bins_s);
m1 = mean(Qpost1/bins_s);
m2 = mean(Qpost2/bins_s);
m3 = mean(Qpost3/bins_s);

OUT.FRates_base = mbase;
OUT.FRates_post1 = m1;
OUT.FRates_post2 = m2;
OUT.FRates_post3 = m3;

OUT.Qbase = Qbase;
OUT.Qbase_bins_uS = bins_uS(:,1) + binsize_Q_ms*1000/2;
OUT.Qpost1 = Qpost1;
OUT.Qpost1_bins_uS = x1(:,1) + binsize_Q_ms*1000/2;
OUT.Qpost2 = Qpost2;
OUT.Qpost2_bins_uS = x2(:,1) + binsize_Q_ms*1000/2;
OUT.Qpost3 = Qpost3;
OUT.Qpost3_bins_uS = x3(:,1) + binsize_Q_ms*1000/2;
OUT.Qall = Qall;
OUT.Qall_x_uS = Qall_x_uS(:,1) + binsize_Q_ms*1000/2;

[OUT.R_base,OUT.R_p_base] = corr(double(Qbase));
[OUT.R_post1,OUT.R_p_post1] = corr(double(Qpost1));
[OUT.R_post2,OUT.R_p_post2] = corr(double(Qpost2));
[OUT.R_post3,OUT.R_p_post3] = corr(double(Qpost3));

if PLOT_IT
    
    figure
    subplot(5,1,1)
    imagesc(start_end_times(:,1)/60e6,[],Q')
    plot_markers(TIMES.EventStartEndUsec/60e6);
    colorbar
    xlabel('min')
    ylabel('neuron')
    title(sprintf('Binsize %d s',binsize_Q_ms/1000))
    plot_markers_simple(TIMES.EventStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    subplot(5,1,2)
    % same but subtract mean rate at beginning.
    
    imagesc(start_end_times(:,1)/60e6,[],Qn')
    plot_markers_simple(TIMES.EventStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    colorbar
    xlabel('min')
    ylabel('neuron')
    title('Subtract rate before injection.')
    
    subplot(5,1,3)
    plot_confidence_intervals(start_end_times(:,1)/60e6,Q')
    ylabel('Mean Neural activity')
    yyaxis right
    plot_confidence_intervals(start_end_times(:,1)/60e6,Qn',[],[.8 .1 .1])
    plot_markers_simple(TIMES.EventStartEndUsec/60e6);
    ylabel('norm neural activity')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot mean firing rates before and after ketamine.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(5,1,4)
    
    %     plot([mbase;m1;m2]')
    plot(start_end_times(:,1)/60e6,sc(:,1:2));
    axis tight
    %     ylabel('Hz')
    %     legend('pre','p1','p2','FontSize',8)
    %     title('Firing rates')
    subplot(5,1,5)
    h = bar(([m1;m2]- mbase)');
    h(1).FaceColor = 'k';
    h(2).FaceColor = 'b';
    set(gca,'Colormap',jet())
    legend('p1','p2','FontSize',8)
    legend boxoff
    title('rate after subtracting baseline')
    xlabel('Neuron')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot overall difference score in firing rates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    histogram_cowen({[m1 - mbase], [m2 - mbase]},.1);
    legend('p1','p2','FontSize',8)
    legend boxoff
    plot_vert_line_at_zero
    xlabel('Change in Firing Rate')
    [~,p1] = ttest([m1 - mbase]);
    [~,p2] = ttest([m2 - mbase]);
    title(sprintf('ttest: p1 %0.4f, p2 %0.4f',p1,p2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at autocorrelations before and after and the difference between
% do it for the pre-injection period but also for the 2-20 minute and 30-50
% minute post periods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[AC,x] = AutoCorrArray_Abhi(TS,bin_size_ac_ms*1000,acorr_width_ms*1000,[TIMES.BaselineUsec;TIMES.PostInjectionUsec]);
% ACd = AC -  AC(:,:,1);
% cm = lines(5);
OUT.AC_base_pre_early_late = AC;
OUT.AC_x_ms = x/1000;
% OUT.AC_conditions = {'pre' 'post 1' 'post 2'};
% disp('ACorr')
% if PLOT_IT
%     titls = OUT.AC_conditions;
%     CIX = x/1000 > 120;
%     
%     figure
%     for ii = 1:length(titls)
%         subplot(1,length(titls),ii)
%         a = squeeze(AC(:,:,ii));
%         a = a - mean(a(:,CIX),2);
%         imagesc(x/1000,[],a)
%         title(titls{ii})
%     end
%     equalize_color_axes
%     
%     figure
%     titls = {'post 1' 'post 2'};
%     for ii = 1:length(titls)
%         hh(ii) = subplot(1,length(titls),ii);
%         imagesc(x/1000,[],AC(:,:,ii+1)- AC(:,:,1))
%         title(titls{ii})
%         colorbar
%         colormap(parula);
%     end
%     equalize_color_axes(hh);
%     
%     figure
%     for ii = 1:3
%         a = squeeze(AC(:,:,ii));
%         a = a - mean(a(:,CIX),2);
%         plot_confidence_intervals(x/1000,a,[],cm(ii,:));
%     end
%     legend_color_text(OUT.AC_conditions,cm(1:3,:));
%     
%     figure
%     plot_confidence_intervals(x/1000,squeeze(AC(:,:,1)),[],cm(1,:))
%     plot_confidence_intervals(x/1000,squeeze(AC(:,:,3)),[],cm(2,:))
%     plot_horiz_line_at_zero
%     
%     
%     figure
%     plot_confidence_intervals(x/1000,squeeze(ACd(:,:,2)),[],cm(1,:))
%     plot_confidence_intervals(x/1000,squeeze(ACd(:,:,3)),[],cm(2,:))
%     plot_horiz_line_at_zero
%     xlabel('ms')
%     ylabel('Difference from baseline')
%     title('Autocorr after drug relative to baseline')
%     legend_color_text({'p1' 'p2'},{cm(1,:) cm(2,:)});
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at local variance...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('LV')
TSx{1} = Restrict(TS,TIMES.BaselineUsec);
TSx{2} = Restrict(TS,TIMES.PostInjectionUsec(1,:));
TSx{3} = Restrict(TS,TIMES.PostInjectionUsec(2,:));
TSx{4} = Restrict(TS,TIMES.PostInjectionUsec(3,:));

SSS = [];
OUT.LocVar_base_pre_early_late = [];
for ii = 1:length(TSx{1})
    for jj = 1:length(TSx)
        OUT.LocVar_base_pre_early_late(ii,jj) = LocalVariance(diff(TSx{jj}{ii}),10e6);
        SSS{ii,jj} = Spiking_summary_stats(TSx{jj}{ii}/1000);
    end
end
OUT.SSS = SSS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the 20 min period right after LDO inj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(intervals_around_evt_min(:,1)) > 3
    TIMES.PostLDOUsec(1,1) = TIMES.Event1StartEndUsec(2) + intervals_around_evt_min(4,1)*60*1e6;
    TIMES.PostLDOUsec(1,2) = TIMES.Event1StartEndUsec(2) + intervals_around_evt_min(4,2)*60*1e6;
    
    % POS and IMU
    ix5 = POS.Time_uS > TIMES.PostLDOUsec(1) & POS.Time_uS < TIMES.PostLDOUsec(2);
    OUT.PostLDO.Speed = nanmean(POS.speed(ix5));
   
    IX5 = IMU.t_uS > TIMES.PostLDOUsec(1) & IMU.t_uS < TIMES.PostLDOUsec(2);
    
    % Store information for analysis
    OUT.PostLDO.Jerk = nanmean(IMU.absjerk(IX5)) ;
    OUT.PostLDO.JerkPC1 = nanmean(IMU.absjerkpc(IX5));
    OUT.PostLDO.IMU_Speed = nanmean(IMU.speed(IX5));

    % Q matrix 
    edges = TIMES.PostLDOUsec(1,1):binsize_Q_ms*1000:TIMES.PostLDOUsec(1,2);
    [QPostLDO,bins_PL] = Bin_ts_array(TS, edges);
    OUT.PostLDO.Q = QPostLDO;
    OUT.PostLDO.bins_uS = bins_PL(:,1) + binsize_Q_ms*1000/2;
    
    m4 = mean(QPostLDO/bins_s);

    OUT.PostLDO.FRates = m4;
    
    [OUT.PostLDO.R,OUT.PostLDO.R_p] = corr(double(QPostLDO));
    
    % AC
    [AC2,x2] = AutoCorrArray_Abhi(TS,bin_size_ac_ms*1000,acorr_width_ms*1000,[TIMES.PostLDOUsec]);
    OUT.PostLDO.AC = AC2;
    OUT.PostLDO.AC_x_ms = x2/1000;
    
    
    % Local variance
    TSxL = Restrict(TS,TIMES.PostLDOUsec(1,:));

    OUT.PostLDO.LocVar = [];
    OUT.PostLDO.SSS = [];
    for ii = 1:length(TSxL)
        OUT.PostLDO.LocVar(ii) = LocalVariance(diff(TSxL{ii}),10e6);
        OUT.PostLDO.SSS{ii} = Spiking_summary_stats(TSxL{ii}/1000);
    end

else
end

OUT.TM = TIMES;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSDs of each spike train - another way to analyze oscillations -
% the stats are easier to interpret I think.
% The problem is that this takes a long time to run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if SPIKE_PSD
%     fqs = 1:.5:90;
%     OUT.SpikePSD_fqs = fqs;
%     OUT.SpikePSD_Before_Early_Late = [];
%     OUT.SpikePSD_Before_Early_LateZ = [];
%     for iS = 1:length(TS)
%         for jj = 1:length(TSx)
%             %         OUT.SpikePSD_Before_Early_Late(iS,:,jj) = Spike_psd(TSx{jj}{iS}/1000,5,fqs,20);
%             %         OUT.SpikePSD_Before_Early_Late(iS,:,jj) = Spike_psd(TSx{jj}{iS}/1000,5,fqs,[]);
%             [OUT.SpikePSD_Before_Early_Late(iS,:,jj),~,~,OUT.SpikePSD_Before_Early_LateZ(iS,:,jj)] = Spike_psd(TSx{jj}{iS}/1000,5,fqs,25);
%         end
%     end
%     if PLOT_IT
%         figure
%         for ii = 1:3
%             subplot(1,3,ii)
%             imagesc(fqs,[],squeeze(OUT.SpikePSD_Before_Early_Late(:,:,ii)))
%         end
%         equalize_color_axes
%         colorbar_label
%         figure
%         for ii = 1:3
%             subplot(1,3,ii)
%             imagesc(fqs,[],squeeze(OUT.SpikePSD_Before_Early_LateZ(:,:,ii)))
%             title('Z')
%         end
%         equalize_color_axes
%         colorbar_label
%     end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Cross Corrs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% small_bin_ms = 2;
% [Combos,CC,CCsh,xdim,stats, SIGNIF] = CrossCorrShuffle_AllCombos(TS,100,small_bin_ms);
% [Combos,CCb,CCshb,xdim,statsb, SIGNIFb] = CrossCorrShuffle_AllCombos(Restrict(TS,TIMES.BaselineUsec),100,small_bin_ms);
% [Combos,CCp,CCshp,xdim,statsp, SIGNIFp] = CrossCorrShuffle_AllCombos(Restrict(TS,TIMES.PostInjectionUsec),100,small_bin_ms);
% 
% OUT.CC = CC;
% OUT.CCb = CCb;
% OUT.CCp = CCp;
% OUT.CCsh = CCsh;
% OUT.CCshb = CCshb;
% OUT.CCshp = CCshp;
% OUT.CCxdim = xdim;
% OUT.Combos = Combos;
% OUT.CC_signif = SIGNIF;
% OUT.CC_signif_b = SIGNIFb;
% OUT.CC_signif_p = SIGNIFp;
% OUT.CC_small_bin_ms = small_bin_ms;
% 
% if PLOT_IT
%     gix = find(SIGNIF(:,2) == 1);
%     for ii = 1:length(gix)
%         iC = gix(ii);
%         figure
%         subplot(1,2,1)
%         stairs(xdim-small_bin_ms/2,CC(iC,:),'LineWidth',3)
%         hold on
%         stairs(xdim-small_bin_ms/2,CCsh(iC,:),'LineWidth',3,'Color','k')
%         xlabel('ms')
%         axis tight
%         pubify_figure_axis
%         plot_vert_line_at_zero
%         ylabel('count')
%         title(sprintf('mono %d glob p %1.3f',stats{iC}.monosynaptic_is_above,stats{iC}.global_stats(1)))
%         
%         subplot(1,2,2)
%         stairs(xdim-small_bin_ms/2,CCb(iC,:),'LineWidth',3)
%         hold on
%         stairs(xdim-small_bin_ms/2,CCp(iC,:),'LineWidth',3)
%         xlabel('ms')
%         axis tight
%         pubify_figure_axis
%         plot_vert_line_at_zero
%         ylabel('count')
%         if ~isempty(statsp{iC})
%             title(sprintf('mono %d glob p %1.3f',statsp{iC}.monosynaptic_is_above,statsp{iC}.global_stats(1)))
%         end
%     end
% end
ca
% [AC,x] = AutoCorrArray(TS,bin_size_ac_ms*1000,acorr_width_ms*1000,[TIMES.BaselineUsec;TIMES.PostInjectionUsec]);
% ACd = AC -  AC(:,:,1);
%  cm = lines(5);
% OUT.AC = AC;
% OUT.AC_x_ms = x/1000;
% OUT.AC_conditions = {'pre' 'post 1' 'post 2'};