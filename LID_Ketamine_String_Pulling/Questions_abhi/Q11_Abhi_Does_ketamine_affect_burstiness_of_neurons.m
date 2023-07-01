function OUT = Q11_Abhi_Does_ketamine_affect_burstiness_of_neurons(type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    type = 'LDOPA';
end
OUT = [];
bin_size_ac_ms = 2;
% binsize_small_ms = 100;
binsize_Q_ms = 100;
acorr_width_ms = 200;    
PLOT_IT = false;

% Load important things
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;
OUT.SP = SP;
OUT.binsize_Q_ms = binsize_Q_ms;
OUT.bin_size_ac_ms = bin_size_ac_ms;

OUT.DEPTHS = DEPTHS; % Send this info out of the function for meta analysis.
OUT.META = META;
OUT.EVT = EVT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the relevant times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_start_end_uS = [];
switch type
    case 'ketamine'
        e_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
        e_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 5 25; 40 60];
        big_peri_event_min = [-25 60];
    case 'LDOPA'
        e_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
        e_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 35 55; 65 85];
        big_peri_event_min = [-25 85];
        
end
if any(E.EventID == 'KetInjectionStart') && any(E.EventID == 'LDOPAInjectionStart')
    disp('Diff ket to ldopa:')
    disp(E.MinFromStart(E.EventID == 'KetInjectionStart') -  E.MinFromStart(E.EventID == 'LDOPAInjectionStart')) 
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

TIMES.PeriKetamineUsec = [TIMES.EventStartEndUsec(1) + big_peri_event_min(1)*60*1e6 TIMES.EventStartEndUsec(1) + big_peri_event_min(2)*60*1e6 ];

OUT.TM = TIMES;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overview: First look at firing rate differences before burstiness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Q,start_end_times] = Bin_ts_array(TS, binsize_Q_ms*1000); % workhorse function for neuron x time matrix
IX = start_end_times(:,1)<TIMES.EventStartEndUsec(1);
Qn = Q-mean(Q(IX,:));
[~,sc] = pca(Qn);

edges = TIMES.BaselineUsec(1,1):binsize_Q_ms*1000:TIMES.BaselineUsec(1,2);
[Qbase,bins_uS] = Bin_ts_array(TS, edges);
edges = TIMES.PostInjectionUsec(1,1):binsize_Q_ms*1000:TIMES.PostInjectionUsec(1,2);
[Qpost1,x1] = Bin_ts_array(TS, edges);
edges = TIMES.PostInjectionUsec(2,1):binsize_Q_ms*1000:TIMES.PostInjectionUsec(2,2);
[Qpost2,x2] = Bin_ts_array(TS, edges);

edges = TIMES.PeriKetamineUsec(1):binsize_Q_ms*1000:TIMES.PeriKetamineUsec(2);
[Qall,Qall_x_uS] = Bin_ts_array(TS, edges);
Qall = int16(Qall);
TSaligned = cell(size(TS));
for ii = 1:length(TS)
    if ~isempty(TS{ii})
        TSaligned{ii} = TS{ii} - TIMES.EventStartEndUsec(1);
    end
end
OUT.TSaligned = TSaligned;

bins_s = binsize_Q_ms/1000;
mbase = mean(Qbase/bins_s);
m1 = mean(Qpost1/bins_s);
m2 = mean(Qpost2/bins_s);

OUT.FRates_base = mbase';
OUT.FRates_post1 = m1';
OUT.FRates_post2 = m2';

OUT.Qbase = Qbase;
OUT.Qbase_bins_uS = bins_uS(:,1) + binsize_Q_ms*1000/2;
OUT.Qpost1 = Qpost1;
OUT.Qpost1_bins_uS = x1(:,1) + binsize_Q_ms*1000/2;
OUT.Qpost2 = Qpost2;
OUT.Qpost2_bins_uS = x2(:,1) + binsize_Q_ms*1000/2;
OUT.Qall = Qall;
OUT.Qall_x_uS = Qall_x_uS(:,1) + binsize_Q_ms*1000/2;

[OUT.R_base,OUT.R_p_base] = corr(double(Qbase));
[OUT.R_post1,OUT.R_p_post1] = corr(double(Qpost1));
[OUT.R_post2,OUT.R_p_post2] = corr(double(Qpost2));

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

[AC,x] = AutoCorrArray(TS,bin_size_ac_ms*1000,acorr_width_ms*1000,[TIMES.BaselineUsec;TIMES.PostInjectionUsec]);
ACd = AC -  AC(:,:,1);
cm = lines(5);
OUT.AC = AC;
OUT.AC_x_ms = x/1000;
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
    title('Autocorr after drug relative to baseline')
    legend_color_text({'p1' 'p2'},{cm(1,:) cm(2,:)});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at local variance, i.e. burstiness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('LV')
TSx{1} = Restrict(TS,TIMES.BaselineUsec);
TSx{2} = Restrict(TS,TIMES.PostInjectionUsec(1,:));
TSx{3} = Restrict(TS,TIMES.PostInjectionUsec(2,:));

SSS = [];
OUT.LocVar_Before_Early_Late = [];
for ii = 1:length(TSx{1})
    for jj = 1:length(TSx)
        OUT.LocVar_Before_Early_Late(ii,jj) = LocalVariance(diff(TSx{jj}{ii}),10e6);
        SSS{ii,jj} = Spiking_summary_stats(TSx{jj}{ii}/1000);
    end
end
OUT.SSS = SSS;

% figure
% histogram_cowen({OUT.LocVar_Before_Early_Late(10,1),OUT.LocVar_Before_Early_Late(10,2),OUT.LocVar_Before_Early_Late(10,3)})
