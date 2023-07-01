function OUT = Q2_What_does_ketamine_do_ensemble()
% Run this function in the data directory.
% This is an example of many different analyses that could be performed but
% it is far too many for a single function. It just demonstrates the
% mixture of things that can be done.
% 
% If you do this for real, focus the function to address just one or a few
% hypotheses. If not, then the code gets too hard to understand and debug.
% 
% Determine if single-units oscillate more after ketamine injection.
% Compare autocorrs before and after injection.
% Cowen 2020
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = [];
bin_size_ms = 2;
acorr_width_ms = 200;       % defining bin sizes and things
PLOT_IT = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function loads a ton of things that are typically used
% for any analysis. It simplifies things.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;
OUT.DEPTHS = DEPTHS;
OUT.META = META;
OUT.EVT = EVT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the relevant times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TIMES.KetStartEndUsec(1,1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
TIMES.KetStartEndUsec(1,2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
TIMES.BaselineUsec(1,1) = TIMES.KetStartEndUsec(1) - 20*60*1e6;
TIMES.BaselineUsec(1,2) = TIMES.KetStartEndUsec(1) - 5*60*1e6;
TIMES.PostInjectionUsec(1,1) = TIMES.KetStartEndUsec(2) + 5*60*1e6;
TIMES.PostInjectionUsec(1,2) = TIMES.KetStartEndUsec(2) + 20*60*1e6;
TIMES.PostInjectionUsec(2,1) = TIMES.KetStartEndUsec(2) + 25*60*1e6;
TIMES.PostInjectionUsec(2,2) = TIMES.KetStartEndUsec(2) + 40*60*1e6;

OUT.TM = TIMES;
%% Load inertial, POS, and calculate string pull rate.
IMU = LK_Load_and_Process_IMU('Inertial_data.mat');
POS = LK_Load_and_Clean_POS('POS.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overview: Look at general movement to determine if ketamine induced motor activity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LK_plot_movement_before_after(
figure
plot(IMU.t_uS/60e6,IMU.absjerk)
hold on
yyaxis right
plot(POS.Time_uS/60e6,POS.xy(:,1))
hold on
plot(POS.Time_uS/60e6,POS.xy(:,2))

IX1 = POS.Time_uS > TIMES.BaselineUsec(1) & POS.Time_uS < TIMES.BaselineUsec(2);
IX2 = POS.Time_uS > TIMES.PostInjectionUsec(1,1) & POS.Time_uS < TIMES.PostInjectionUsec(1,2);
IX3 = POS.Time_uS > TIMES.PostInjectionUsec(2,1) & POS.Time_uS < TIMES.PostInjectionUsec(2,2);
OUT.Speed_Before_Early_Late = [nanmean(POS.speed(IX1)) nanmean(POS.speed(IX2)) nanmean(POS.speed(IX3))];

IX1 = IMU.t_uS > TIMES.BaselineUsec(1) & IMU.t_uS < TIMES.BaselineUsec(2);
IX2 = IMU.t_uS > TIMES.PostInjectionUsec(1,1) & IMU.t_uS < TIMES.PostInjectionUsec(1,2);
IX3 = IMU.t_uS > TIMES.PostInjectionUsec(2,1) & IMU.t_uS < TIMES.PostInjectionUsec(2,2);
OUT.Jerk_Before_Early_Late = [nanmean(IMU.absjerk(IX1)) nanmean(IMU.absjerk(IX2)) nanmean(IMU.absjerk(IX3))];
OUT.JerkPC1_Before_Early_Late = [nanmean(IMU.absjerkpc(IX1)) nanmean(IMU.absjerkpc(IX2)) nanmean(IMU.absjerkpc(IX3))];



% V = [POS.speed(IX1); POS.speed(IX2); POS.speed(IX3)];
% G = [ones(sum(IX1),1);2*ones(sum(IX2),1);3*ones(sum(IX3),1)];



figure
subplot(3,1,1)
bar(OUT.Speed_Before_Early_Late)
set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
ylabel('Speed')
pubify_figure_axis
subplot(3,1,2)
bar(OUT.Jerk_Before_Early_Late)
set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
ylabel('Jerk')
pubify_figure_axis
subplot(3,1,3)
bar(OUT.JerkPC1_Before_Early_Late)
set(gca,'XTickLabel',{'before' 'after 1' 'after 2'})
ylabel('JerkPC')
pubify_figure_axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overview: Generate a neuron by time matrix and mark when ketamine was delivered.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binsize_ms = 10000;
[Q,start_end_times] = Bin_ts_array(TS, binsize_ms*1000);
IX = start_end_times(:,1)<TIMES.KetStartEndUsec(1);
Qn = Q-mean(Q(IX,:));
[Qbase,bins_uS] = Bin_ts_array(Restrict(TS,TIMES.BaselineUsec), binsize_ms*1000);
[Qpost1,x1] = Bin_ts_array(Restrict(TS,TIMES.PostInjectionUsec(1,:)), binsize_ms*1000);
[Qpost2,x2] = Bin_ts_array(Restrict(TS,TIMES.PostInjectionUsec(2,:)), binsize_ms*1000);
bins_s = binsize_ms/1000;
mbase = mean(Qbase/bins_s);
m1 = mean(Qpost1/bins_s);
m2 = mean(Qpost2/bins_s);

OUT.FRates_base = mbase;
OUT.FRates_post1 = m1;
OUT.FRates_post2 = m2;

OUT.Qbase = Qbase;
OUT.Qbase_bins_uS = bins_uS;
OUT.Qpost1 = Qpost1;
OUT.Qpost1_bins_uS = x1;
OUT.Qpost2 = Qpost2;
OUT.Qpost2_bins_uS = x2;

if PLOT_IT
    
    figure
    subplot(5,1,1)
    imagesc(start_end_times(:,1)/60e6,[],Q')
    plot_markers(TIMES.KetStartEndUsec/60e6);
    colorbar
    xlabel('min')
    ylabel('neuron')
    title(sprintf('Binsize %d s',binsize_ms/1000))
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    subplot(5,1,2)
    % same but subtract mean rate at beginning.
    
    imagesc(start_end_times(:,1)/60e6,[],Qn')
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    plot_markers_simple(TIMES.BaselineUsec/60e6,[],[],'w');
    plot_markers_simple(TIMES.PostInjectionUsec/60e6,[],[],'g');
    
    colorbar
    xlabel('min')
    ylabel('neuron')
    title('Subtract rate before injection.')
    
    subplot(5,1,3)
    plot_confidence_intervals(start_end_times(:,1)/60e6,Q')
    yyaxis right
    plot_confidence_intervals(start_end_times(:,1)/60e6,Qn',[],[.8 .1 .1])
    plot_markers_simple(TIMES.KetStartEndUsec/60e6);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot mean firing rates before and after ketamine.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(5,1,4)
    
    bar([mbase;m1;m2]')
    ylabel('Hz')
    legend('pre','p1','p2','FontSize',8)
    title('Firing rates')
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

[AC,x] = AutoCorrArray(TS,bin_size_ms*1000,acorr_width_ms*1000,[TIMES.BaselineUsec;TIMES.PostInjectionUsec]);
ACd = AC -  AC(:,:,1);
cm = lines(5);
OUT.AC = AC;
OUT.AC_x_ms = x/1000;
OUT.AC_conditions = {'pre' 'post 1' 'post 2'};
disp('ACorr')
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
    title('Autocorr after ketamine relative to baseline')
    legend_color_text({'p1' 'p2'},{cm(1,:) cm(2,:)});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smooth autocorr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SACorr')
TSx{1} = Restrict(TS,TIMES.BaselineUsec);
TSx{2} = Restrict(TS,TIMES.PostInjectionUsec(1,:));
TSx{3} = Restrict(TS,TIMES.PostInjectionUsec(2,:));
binsize_ms = 2;
smth_win_ms = 8;
acorr_size_ms = 250;
jitter_ms = 500;
for ii = 1:length(TS)
    for jj = 1:length(TSx)
         [OUT.AcorrSmth_Before_Early_Late(ii,:,jj),OUT.AcorrSmth_x_ms, CI] = ...
             CrossCorr_smooth(TSx{jj}{ii}/1000,[],binsize_ms,acorr_size_ms/binsize_ms,smth_win_ms/binsize_ms,jitter_ms);
         OUT.AcorrSmthCI95up_Before_Early_Late(ii,:,jj)= CI(1,:);
         OUT.AcorrSmthCI95down_Before_Early_Late(ii,:,jj)= CI(2,:);
         %         CrossCorr_smooth(TSx{jj}{ii}/1000,[],binsize_ms,acorr_size_ms/binsize_ms,smth_win_ms/binsize_ms,jitter_ms)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at local variance...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('LV')
SSS = [];
for ii = 1:length(TS)
    for jj = 1:length(TSx)
        OUT.LocVar_Before_Early_Late(ii,jj) = LocalVariance(diff(TSx{jj}{ii}),10e6);
        SSS = Spiking_summary_stats(TSx{jj}{ii}/1000);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at Cross Corrs between neurons. Only store ones that are significnat.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSDs of each spike train - another way to analyze oscillations -
% the stats are easier to interpret I think.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    fqs = 1:.5:70;
    OUT.SpikePSD_fqs = fqs;
    OUT.SpikePSD_Before_Early_Late = [];
    OUT.SpikePSD_Before_Early_LateZ = [];
    for iS = 1:length(TS)
        for jj = 1:length(TSx)
            %         OUT.SpikePSD_Before_Early_Late(iS,:,jj) = Spike_psd(TSx{jj}{iS}/1000,5,fqs,20);
            %         OUT.SpikePSD_Before_Early_Late(iS,:,jj) = Spike_psd(TSx{jj}{iS}/1000,5,fqs,[]);
            [OUT.SpikePSD_Before_Early_Late(iS,:,jj),~,~,OUT.SpikePSD_Before_Early_LateZ(iS,:,jj)] = Spike_psd(TSx{jj}{iS}/1000,5,fqs,25);
        end
    end
    figure
    for ii = 1:3
        subplot(1,3,ii)
        imagesc(fqs,[],squeeze(OUT.SpikePSD_Before_Early_Late(:,:,ii)))
    end
    equalize_color_axes
    colorbar_label
    figure
    for ii = 1:3
        subplot(1,3,ii)
        imagesc(fqs,[],squeeze(OUT.SpikePSD_Before_Early_LateZ(:,:,ii)))
        title('Z')
    end
    equalize_color_axes
    colorbar_label
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at higher-order measures of synchrony BETWEEN neurons like
% correlations, coincidences (within 20ms), entropy.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a more finely binned spiketime matrix and calculate the degree of
% correlation bettween neurons...
binsize_ms = 25;
window_size_sec = 120;
window_size_bins = window_size_sec/(binsize_ms/1000);

[Q,Q_x_uS] = Bin_ts_array(TS,binsize_ms*1000);

[self_similarity, ss_uS ]= Sliding_r(Q, Q_x_uS(1:end,1), window_size_bins);

IX1 = ss_uS(:,1) > TIMES.BaselineUsec(1) & ss_uS(:,1) < TIMES.BaselineUsec(2);
IX2 = ss_uS(:,1) > TIMES.PostInjectionUsec(1,1) & ss_uS(:,1) < TIMES.PostInjectionUsec(1,2);
IX3 = ss_uS(:,1) > TIMES.PostInjectionUsec(2,1) & ss_uS(:,1) < TIMES.PostInjectionUsec(2,2);
OUT.SelfSim_Before_Early_Late = [nanmean(self_similarity(IX1)) nanmean(self_similarity(IX2)) nanmean(self_similarity(IX3))];

OUT.self_sim = self_similarity;
OUT.self_sim_x_sec = ss_uS;

if PLOT_IT
    figure
    plot(ss_uS(:,1)/60e6,self_similarity)
    plot_vert_line_at_zero(TIMES.KetStartEndUsec/1e6/60);
    axis tight
    title('Self-similarity')
    xlabel('min')
end

%%%%%%%% ENTROPY: RANDOMNESS IN THE POPULATION CODE. %%%%%%%%%%%%
C = kmeans(Q, 10);
[ent,time_ts]= Sliding_entropy([Q_x_uS(:) C(:)], window_size_bins, window_size_bins); % 10*20,5*20
x_min = linspace(5,META.sec_of_recording-5,length(ent))/60;

IX1 = x_min*60e6 > TIMES.BaselineUsec(1) & x_min*60e6 < TIMES.BaselineUsec(2);
IX2 = x_min*60e6 > TIMES.PostInjectionUsec(1,1) & x_min*60e6 < TIMES.PostInjectionUsec(1,2);
IX3 = x_min*60e6 > TIMES.PostInjectionUsec(2,1) & x_min*60e6 < TIMES.PostInjectionUsec(2,2);
OUT.Entropy_Before_Early_Late = [nanmean(ent(IX1)) nanmean(ent(IX2)) nanmean(ent(IX3))];


OUT.entropy = ent;
OUT.entropy_x_min = x_min;

if PLOT_IT
    figure
    plot(x_min,ent)
    plot_vert_line_at_zero(TIMES.KetStartEndUsec/1e6/60);
    axis tight
    title('Entropy')
    xlabel('min')
end

%%%%%% Coincidences....
binsize_ms = 20;
bins_uS = 0:binsize_ms*1e3:META.sec_of_recording*1e6;
[Q] = histcounts_cowen(TS,bins_uS);
nc = sum(Q);
ncd = decimate(nc,100);
ncdgz = decimate(double(nc>0),100);
x_uS = linspace(bins_uS(1),bins_uS(end),length(ncd));

IX1 = x_uS > TIMES.BaselineUsec(1) & x_uS < TIMES.BaselineUsec(2);
IX2 = x_uS > TIMES.PostInjectionUsec(1,1) & x_uS < TIMES.PostInjectionUsec(1,2);
IX3 = x_uS > TIMES.PostInjectionUsec(2,1) & x_uS < TIMES.PostInjectionUsec(2,2);
OUT.Coincidences_Before_Early_Late = [nanmean(ncd(IX1)) nanmean(ncd(IX2)) nanmean(ncd(IX3))];
OUT.n_coincidences = ncd;
OUT.n_coincidences_x_uS = x_uS;

%%%%%% Network Complexity
binsize_ms = 20;
[Q,Q_x_uS] = Bin_ts_array(TS,binsize_ms*1000);
Q_x_uS = mean(Q_x_uS,2);
% [ND]= Sliding_dynamics(Q,20000/binsize_ms);
%
IX1 = Q_x_uS > TIMES.BaselineUsec(1) & Q_x_uS < TIMES.BaselineUsec(2);
IX2 = Q_x_uS > TIMES.PostInjectionUsec(1,1) & Q_x_uS < TIMES.PostInjectionUsec(1,2);
IX3 = Q_x_uS > TIMES.PostInjectionUsec(2,1) & Q_x_uS < TIMES.PostInjectionUsec(2,2);

OUT.nEffDim_Before_Early_Late(1) = n_effective_dimensions(Q(IX1,:));
OUT.nEffDim_Before_Early_Late(2) = n_effective_dimensions(Q(IX2,:));
OUT.nEffDim_Before_Early_Late(3) = n_effective_dimensions(Q(IX3,:));



if PLOT_IT
    figure
    plot(x_uS/60e6, conv(ncd,hanning(100),'same'));
    axis tight
    plot_vert_line_at_zero(TIMES.KetStartEndUsec/1e6/60);
    title('nCoincidences')
    
    figure
    plot(x_uS/60e6, conv(ncdgz,hanning(100),'same'));
    axis tight
    title('n more than 0 Coincidence')
    plot_vert_line_at_zero(TIMES.KetStartEndUsec/1e6/60);
end
%%%%% junk
if 0
    
    binsize_ms = 10000;
    bins = 0:binsize_ms/1000:META.sec_of_recording;
    window_size_sec = 120;
    window_size_bins = window_size_sec/(binsize_ms/1000);
    Q = nan(length(SP),length(bins)-1);
    for iS = 1:length(SP)
        Q(iS,:) = histcounts(SP(iS).t_uS/1e6,bins);
    end
    
    
    
    %%
    [CC,x] = CrossCorr(SP(1).t_uS/1000,SP(2).t_uS/1000,bin_size_ms,acorr_width_ms/bin_size_ms);
    figure
    plot(x,CC)
    xlabel('ms');
    ylabel('Hz')
    
    for ii = 1:length(SP)
        Tbefore_s = Restrict(SP(ii).t_uS/1e6,[0 TIMES.KetStartEndUsec(1)/1e6 - 20]); % 20 seconds before ket injection
        Tafter_s = Restrict(SP(ii).t_uS/1e6,[TIMES.KetStartEndUsec(2)/1e6 + 60, 180*60]); % 60 seconds after ket injection
        [ACbef(ii,:),x_ms] = AutoCorr(Tbefore_s*1000,bin_size_ms,acorr_width_ms/bin_size_ms);
        [ACaft(ii,:),x_ms] = AutoCorr(Tafter_s*1000,bin_size_ms,acorr_width_ms/bin_size_ms);
    end
    
    figure
    subplot(2,2,1)
    imagesc(x_ms,[],ACbef)
    title('Before')
    subplot(2,2,2)
    imagesc(x_ms,[],ACaft)
    title('After')
    equalize_color_axes
    subplot(2,2,3:4)
    plot_confidence_intervals(x_ms,ACbef,[],GP.Colors.Control)
    plot_confidence_intervals(x_ms,ACaft,[],GP.Colors.Ketamine)
    xlabel('ms')
end