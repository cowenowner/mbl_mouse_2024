function OUT = Q3_does_stim_affect_ensemble_act()
%%
%
% Assumes that the Neuropixels data has been synchronized with Tprime
% (see NPXL_... functions) and an event file has been created.
%
% This does  NOT analyze the DA data - for example when no DA was recorded.
%
% Assumes the following files are in the dataset directory:
% AllSpikes.mat
% Events.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cowen 2023 (may not work on the 2022 data - change of format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
binsize_ms = 500; % note - for analyses, we should really ignore bins that overlap the scans since those bins will always be zeros.
PLOT_IT = false;
selfsim_win_size_bins = 20;
plot_each_peth = 0;
time_before_sec = 25;
time_after_sec = 60*2;
peth_bin_ms = 550;
clrs = lines(200);

[~,set_name] = fileparts(pwd);

[p,rat_ID] = fileparts(fileparts(pwd));
rat_num = str2double(rat_ID(4:end));
% OUT.rat = num2str(rat_id(4:end))

load('Events.mat','EVT')
% DANA_visualize_Events(EVT)

switch rat_ID
    case 'Rat000'
        all_scan_times = EVT.scan_time.time_sec; % works for rat 425.
        stim_times_sec = EVT.stim_time.time_sec;
    case 'Rat425'
        all_scan_times = EVT.sync3.time_sec; % works for rat 425.
        stim_times_sec = EVT.stim.time_sec;
    case {'Rat438' 'Rat439' 'Rat445' }
        all_scan_times = EVT.fscv.time_sec; % a bug - this is the case for rat 438 439
        stim_times_sec = EVT.stim.time_sec;
end

% block_se = [EVT.fscv.time_sec(:) EVT.fscv.time_sec(:) + block_dur_sec];
[SP] = DANA_process_spikes();
T_uS = {SP.t_uS}; WV = [SP.WaveformBest]';

% figure;imagesc(sort_matrix(WV,'pc1'))
% Verify...
% figure;plot_raster(T_uS,[],'time_divisor',60e6);xlabel('minutes')
% figure;Plot_spike_by_time_matrix(T_uS, .1e6, 'time_divisor',1e6,'xlab','sec');
%%%%%%%%%%%%%%%%%
% Find valid scan times (the scan times that only occurred within a
% recording block).
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
% Divine trial information from the stim times.
% Determine the stim block start and end times.
%%%%%%%%%%%%%%%%%
TI = DANA_trial_info_from_stim_times(stim_times_sec);
% TI.Hz(10:end);
OUT.TI = TI;
OUT.PETH_x_axis_ms = [];
%% Compute StimTAs aligned to stim offset for the different conditions.

% Create table of aligned responses.
[OUT.TBL, OUT.INFO] = DANA_Condition_Aligned_Spike_Table_From_TI(SP,TI,rat_num,'time_after_sec',time_after_sec);

% Process the ensemble data: Neuron x time matrices.
if 0
    mn = mean(ALL_M(:,BEFIX),2);
    sd = std(ALL_M(:,BEFIX),[],2);
    ALL_M_z = (ALL_M-mn)./sd;

    mn = mean(ALL_Mnorm(:,BEFIX),2);
    sd = std(ALL_Mnorm(:,BEFIX),[],2);
    ALL_M_norm_z = (ALL_Mnorm-mn)./sd;
end
[Q, edges_uS] = histcounts_cowen(T_uS, 'binsize', binsize_ms*1000);
Qz = Z_scores(Q);
Q_uS = edges_uS(1:end-1) + (edges_uS(2) - edges_uS(1))/2;

% OUT.Q = Q;
% OUT.Q_uS = Q_uS;
% % OUT.ALL_M = ALL_M;
% % OUT.ALL_Mnorm = ALL_Mnorm;
% OUT.ALL_Hz_G = ALL_Hz_G;
% OUT.ALL_LV_G = ALL_LV_G;
% OUT.Hz = Hz;
% OUT.LV = LV;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the data.
% big picture plot as an overview....
if PLOT_IT
    figure
    plot_raster(T_uS,[],'time_divisor',3600e6);xlabel('hours')
    xlabel('Hours'); ylabel('Neuron')
    hold on
    % plot(stim_times_sec/3600,ones(size(stim_times_sec)),'w^')
    plot_markers_simple(TI.end_times_sec(TI.Hz_group == 10)/3600,[],1,[1,0,0])
    plot_markers_simple(TI.start_times_sec(TI.Hz_group == 10)/3600,[],1,[1,0,0])
    plot_markers_simple(TI.end_times_sec( TI.Hz_group == 20)/3600,[],1,[0,1,0])
    plot_markers_simple(TI.start_times_sec( TI.Hz_group == 20)/3600,[],1,[0,1,0])
    plot_markers_simple(TI.end_times_sec(TI.Hz_group == 60)/3600,[],1,[1,1,1])
    set(gcf,'Position',[20         115        1392         689])
    text(TI.end_times_sec/3600, Cols(Qz)*ones(1,length(TI.end_times_sec))+1, num2str(TI.LV_group),'FontSize',8)
    text(TI.end_times_sec/3600, Cols(Qz)*ones(1,length(TI.end_times_sec))+2, num2str(TI.Hz_group),'FontSize',8)
    sgtitle(sprintf('%s',rat_ID))
% return

    %
    tmp_clrs = [1 .2 .2;.2 .2 1];
    figure
    uHz = unique(OUT.TBL.Hz);
    for ii = 1:length(uHz)
        IX = OUT.TBL.Hz == uHz(ii);
        plot_confidence_intervals(OUT.INFO.PETH_x_axis_ms/1000,OUT.TBL.mean_PETHnorm(IX,:),[],tmp_clrs(ii,:))
    end
    title('Hz'); ylabel('Norm rate'); xlabel('sec')
    pubify_figure_axis
    legend_color_text({'10' '20'},tmp_clrs(1:2,:))


    figure
    uLVs = unique(OUT.TBL.LV);
    for ii = 1:length(uLVs)
        IX = OUT.TBL.LV == uLVs(ii);
        plot_confidence_intervals(OUT.INFO.PETH_x_axis_ms/1000,OUT.TBL.mean_PETHnorm(IX,:),[],clrs(ii,:))
    end
    title('LV'); ylabel('Norm rate'); xlabel('sec')
    pubify_figure_axis
    legend_color_text({'0' '.3' '1' '1.2'},clrs(1:4,:))

    % Compare correlations before and after.

    figure
    for ii = 1:length(uHz)
        IX = OUT.TBL.Hz == uHz(ii);
        M = OUT.TBL.LocalVarianceAroundStim(IX,:);
        M = M - M(:,1);
        subplot(1,2,ii)
        boxplot(M,'notch','on')
        title(sprintf('%d Hz',uHz(ii)))
    end

    figure
    for ii = 1:length(uHz)
        IX = OUT.TBL.Hz == uHz(ii);
        M = OUT.TBL.CVAroundStim(IX,:);
        M = M - M(:,1);
        subplot(1,2,ii)
        boxplot(M,'notch','on')
        ylabel('CV')
        title(sprintf('%d Hz',uHz(ii)))
    end


    figure
    for ii = 1:length(uLVs)
        IX = OUT.TBL.LV == uLVs(ii);
        M = OUT.TBL.LocalVarianceAroundStim(IX,:);
        M = M - M(:,1);

        subplot(1,length(uLVs),ii)
        boxplot(M,'notch','on')
        ylabel('LV')
        title(sprintf('%1.2f LV',uLVs(ii)))
    end


    figure
    for ii = 1:length(uLVs)
        IX = OUT.TBL.LV == uLVs(ii);
        M = OUT.TBL.CVAroundStim(IX,:);
        M = M - M(:,1);
        subplot(1,length(uLVs),ii)
        boxplot(M,'notch','on')
        ylabel('CV')
        title(sprintf('%1.2f LV',uLVs(ii)))
    end
end

FRATE_NORM = OUT.TBL.mean_PETHnorm;
% FRATE_NORM = ALL_M_z;

if 0
    clrs = lines;
    figure
    subplot(4,2,1:2)
    % imagesc(Q_uS/3600e6,[],Qz')
    plot_raster(T_uS,[],'time_divisor',3600e6);xlabel('hours')
    % clim([-1 3])
    xlabel('Hours'); ylabel('Neuron')
    hold on
    plot(stim_times_sec/3600,ones(size(stim_times_sec)),'w^')
    plot_markers_simple(TI.end_times_sec(TI.Hz_group == 10)/3600,[],1,[1,0,0])
    plot_markers_simple(TI.end_times_sec( TI.Hz_group == 20)/3600,[],1,[0,1,0])
    plot_markers_simple(TI.end_times_sec(TI.Hz_group == 60)/3600,[],1,[1,1,1])

    text(TI.end_times_sec/3600, Cols(Qz)*ones(1,length(TI.end_times_sec)), num2str(TI.LV_group))
    text(TI.end_times_sec/3600, zeros(1,length(TI.end_times_sec)), num2str(TI.Hz_group))

    % yyaxis right
    % plot(cv_data(:,1)/3600,cv_data(:,2),'w.')
    % ylabel('[DA]')
    % title(sprintf('%s Population activity (R = 10Hz, G = 20Hz)',set_name))

    subplot(4,2,3)

    IX = ALL_Hz_G == 10;
    lv = ALL_LV_G(IX);
    m  = ALL_M(IX(:),:);
    [~,six] = sort(lv);
    m = m(six,:);
    imagesc(x_axis_ms/1000,[],m);
    title('10Hz')
    xlabel('s'); ylabel('Trial sorted on LV')
    plot_vert_line_at_zero;
    plot_vert_line_at_zero(-10);
    cax = caxis;


    subplot(4,2,4)
    IX = ALL_Hz_G == 20;
    lv = ALL_LV_G(IX);
    m  = ALL_M(IX(:),:);
    [~,six] = sort(lv);
    m = m(six,:);
    imagesc(x_axis_ms/1000,[],m);
    title('20Hz')
    xlabel('s'); ylabel('Trial')
    plot_vert_line_at_zero;
    plot_vert_line_at_zero(-10);
    caxis(cax);

    subplot(4,2,5:6)
    plot_confidence_intervals(x_axis_ms/1000,FRATE_NORM(ALL_Hz_G == 10,:),[],[1,0,0])
    plot_confidence_intervals(x_axis_ms/1000,FRATE_NORM(ALL_Hz_G == 20,:),[],[0,1,0])
    plot_vert_line_at_zero
    plot_vert_line_at_zero(-10);

    legend_color_text({'10' '20'},{[1,0,0] [0,1,0]})
    xlabel('s')
    ylabel('sd from base')
    title('Rate')

    subplot(4,2,7:8)
    LVs = unique(ALL_LV_G);
    clr = lines(5);
    for ii = 1:length(LVs)
        plot_confidence_intervals(x_axis_ms/1000,FRATE_NORM(ALL_LV_G == LVs(ii),:),[],clrs(ii,:))
    end
    plot_vert_line_at_zero
    plot_vert_line_at_zero(-10);

    legend_color_text({'0' '.3' '1' '1.2'},clrs(1:4,:))
    title('LV')
    xlabel('s')
    ylabel('sd from base')
end

%%
sm_bin = 10;
[pc,sc ] = pca(Q);
[self_sim_r, self_sim_r_x ]= Sliding_r(Q, Q_uS/3600e6, selfsim_win_size_bins);
self_sim_r_x = mean(self_sim_r_x,2);
pop_rate = mean(Q,2);
S = Sparseness_measures(Q);
sted = [Find_start(T_uS) Find_end(T_uS)];
edges = sted(1):10e6:sted(end);
LV = nan(length(T_uS),length(edges)-1);
for iE = 1:(length(edges)-1)
    t = Restrict(T_uS,edges(iE:iE+1));
    for iN = 1:length(t)
        if length(t{iN})>10
            LV(iN,iE) = LocalVariance(diff(t{iN}));
        end
    end
end

% smooth it
pop_rate = movmean(pop_rate,sm_bin,'omitnan');
sc(:,1) = movmean(sc(:,1),sm_bin,'omitnan');
sc(:,2) = movmean(sc(:,2),sm_bin,'omitnan');
S.CV = movmean(S.CV,sm_bin,'omitnan');
S.Prop_Active = movmean(S.Prop_Active,sm_bin,'omitnan');
S.kurtosis = movmean(S.kurtosis,sm_bin,'omitnan');


if PLOT_IT
    figure
    ax = [];
    ax(1) = subplot(4,1,1);
    imagesc(Q_uS/3600e6,[],Qz')
    clim([-1 3])
    xlabel('Hours'); ylabel('Neuron')
    hold on
    plot(stim_times_sec/3600,ones(size(stim_times_sec)),'w^')
    plot_markers_simple(TI.end_times_sec(TI.Hz_group == 10)/3600,[],1,[1,0,0])
    plot_markers_simple(TI.end_times_sec( TI.Hz_group == 20)/3600,[],1,[0,1,0])
    plot_markers_simple(TI.end_times_sec(TI.Hz_group == 60)/3600,[],1,[1,1,1])

    text(TI.end_times_sec/3600, Cols(Qz)*ones(1,length(TI.end_times_sec)), num2str(TI.LV_group))
    text(TI.end_times_sec/3600, zeros(1,length(TI.end_times_sec)), num2str(TI.Hz_group))
    title(sprintf('%s Population activity (R = 10Hz, G = 20Hz)',set_name))

    %
    ax(2) = subplot(4,1,2);

    plot(Q_uS/3600e6,pop_rate)
    ylabel('Pop Rate')
    yyaxis right
    plot(Q_uS/3600e6,sc(:,1))
    hold on
    plot(Q_uS/3600e6,sc(:,2))
    ylabel('PC 1 2')
    axis tight

    ax(3) = subplot(4,1,3);
    plot(edges(1:end-1)/3600e6,mean(LV,'omitnan'))
    % plot(Q_uS/3600e6,S.Prop_Active)
    % ylabel('Prop Active')
    ylabel('Sliding LV')
    yyaxis right
    plot(Q_uS/3600e6,S.kurtosis)
    ylabel('kurtosis')
    axis tight

    ax(4) = subplot(4,1,4);
    plot(self_sim_r_x,self_sim_r)
    ylabel('Adjacent Similarity')
    axis tight

    linkaxes(ax,'x')

    %%
    figure
    subplot(2,2,1:2)
    plot(pc(:,1));hold on;plot(pc(:,2))
    title('PC 1 an 2')
    subplot(2,2,3)
    plot(sc(:,1),sc(:,2),'.')
    subplot(2,2,4)
    plot(sc(:,2),sc(:,3),'.')

    % Sparseness_measures(Q);

    %% Correlations
    window_size_bins = 20;
    [self_similarity, out_times ]= Sliding_r(Q, Q_uS/3600e6, window_size_bins);

    figure
    plot(mean(out_times,2), self_similarity)

end
% % Local variance
% for iN = 1:length(T_uS)
%     LV(iN,:) = LocalVariance(diff(T_uS{iN},))
% end

%% END
