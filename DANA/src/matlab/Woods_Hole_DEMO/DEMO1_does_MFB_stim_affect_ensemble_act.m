% function OUT = Q2_does_DA_affect_ensemble_act_2023()
%%
%
% Assumes that the Neurpixels data has been synchronized with Tprime
% (see NPXL_... functions) and an event file has been created.
% Assumes had-normalized FSCV data in the .csv file.
%
% Assumes the following files are in the dataset directory:
% AllSpikes.mat
% Events.mat
% a folder called WCCV_PCA that has the many .txt tab-separated value files
% and a single .xlsx file that has the order in which the blocks occurred.
% Good sessions:
% GitHub\DANA\Data\Acute\Processed_Data\20221101
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cowen 2023 (may not work on the 2022 data - change of format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd C:\Users\cowen\Documents\GitHub\DANA\Data\Acute\Processed_Data\Rat425\01
% linkaxes([ax1 ax2 ax3],'xy') % this command often causes matlab to crash.

% Works on Rat 445, 439, not on 425.

OUT = [];
binsize_ms = 50; % note - for analyses, we should really ignore bins that overlap the scans since those bins will always be zeros.
sliding_window_ms = 4000; % for the dyanmics analysis and smoothing.
block_dur_sec = 10;
PLOT_IT = true;
if exist('.\WCCV_PCA')
    WCCV_data_folder = '.\WCCV_PCA';
else
    d = dir('Rat_*_i_vs_t');
    WCCV_data_folder = d(1).name;   
end
selfsim_win_size_bins = 20;

[~,set_name] = fileparts(pwd);
[p,rat_ID] = fileparts(fileparts(pwd));
rat_num = str2double(rat_ID(4:end));


load('Events.mat','EVT')
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

block_se = [EVT.fscv.time_sec(:) EVT.fscv.time_sec(:) + block_dur_sec];

[SP] = DANA_process_spikes();
T_uS = {SP.t_uS}; WV = [SP.WaveformBest]';

% Verify...
if 0
    DANA_visualize_Events(EVT)
    Plot_spike_by_time_matrix(T_uS, .1e6);
end
%%%%%%%%%%%%%%%%%
% Find valid scan times (the scan times that only occurred within a
% recording block).
%%%%%%%%%%%%%%%%%
STUPID_WCCV_SHIFT_ASSUMPTION_SEC = 5.01; % so arbitrary, so asynchronous, so wrong.
good_npxl_scan_times_sec = Restrict(all_scan_times,block_se-STUPID_WCCV_SHIFT_ASSUMPTION_SEC);
% block_start_time_sec = EVT.block_time.time_sec(:) - STUPID_WCCV_SHIFT_ASSUMPTION_SEC; % again, stupid that you have to subtract 5. Makes no sense.
% even if we have the start time, how do we know how long each scan is
% since the block end time does not tell you this.

%%%%%%%%%%%%%%%%%
% Divine trial information from the stim times.
% Determine the stim block start and end times.
%%%%%%%%%%%%%%%%%
TI = DANA_trial_info_from_stim_times(stim_times_sec);
NPXL_block_start_sec = TI.start_times_sec - STUPID_WCCV_SHIFT_ASSUMPTION_SEC; % yeah - dumb.

% Load the CV data.
[BLOCK_INFO, CV_DATA ] = DANA_load_wccv_data_from_folder(WCCV_data_folder);
block_se(:,2) - block_se(:,1);
fprintf('Neuropixels recorded %d blocks. WCCV reported %d.\n', Rows(block_se), length(BLOCK_INFO))
fprintf('Neuropixels recorded valid %d scans. WCCV reported %d.\n', Rows(good_npxl_scan_times_sec), Rows(CV_DATA))

[cv_data, BLOCK_INFO, ALLCV ] = DANA_sync_WCCV_with_NPXL_using_TI(BLOCK_INFO, NPXL_block_start_sec, all_scan_times, TI);

% double check to make sure times are aligned.
[BLOCK_INFO(13:end).Hz]
TI.Hz(10:end)

st = Rows(CV_DATA);
CV_DATA(:,5) = good_npxl_scan_times_sec(Rows(CV_DATA));


%% Compute STAs aligned to stim offset for the different conditions.
Hz = [10 20];
LV = [0 .3 1 1.2];
plot_each_peth = 0;
time_before_sec = 25;
time_after_sec = 50;
bin_ms = 550;

ALL_M = [];
ALL_Hz_G = [];
ALL_LV_G = [];
cnt = 1;
for iG = 1:length(Hz)
    for iLV = 1:length(LV)
        ix = TI.Hz_group == Hz(iG) & TI.LV_group == LV(iLV);
        for iN = 1:length(T_uS)
            %         PETH_raster(T_uS{iN}/100,TI.end_times_sec(ix)*10000,20,time_before_sec*1000,time_after_sec*1000);
            [M, x_axis, A_msec, h ] = PETH_raster(T_uS{iN}/100,TI.end_times_sec(ix)*10000,bin_ms,time_before_sec*1000,time_after_sec*1000);
            ALL_M(cnt,:) = mean(M);
            ALL_Hz_G(cnt) = Hz(iG);
            ALL_LV_G(cnt) = LV(iLV);
            cnt = cnt + 1;
        end
    end

end
% Process the ensemble data: Neuron x time matrices.
base_ix = x_axis < -11;
mn = mean(ALL_M(:,base_ix),2);
sd = std(ALL_M(:,base_ix),[],2);
ALL_M_norm = (ALL_M-mn)./sd;

[Q, edges_uS] = histcounts_cowen(T_uS, 'binsize', binsize_ms*1000);
Qz = Z_scores(Q);
Q_uS = edges_uS(1:end-1) + (edges_uS(2) - edges_uS(1))/2;


% Align_and_interpolate_on_events

% Plot the data.
% big picture plot as an overview....
clrs = lines;

figure
subplot(4,2,1:2)
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

yyaxis right
plot(cv_data(:,1)/3600,cv_data(:,2),'w.')
ylabel('[DA]')
title(sprintf('%s Population activity (R = 10Hz, G = 20Hz)',set_name))

subplot(4,2,3)

IX = ALL_Hz_G == 10;
lv = ALL_LV_G(IX);
m  = ALL_M(IX(:),:);
[~,six] = sort(lv);
m = m(six,:);
imagesc(x_axis/1000,[],m);
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
imagesc(x_axis/1000,[],m);
title('20Hz')
xlabel('s'); ylabel('Trial')
plot_vert_line_at_zero;
plot_vert_line_at_zero(-10);
caxis(cax);

subplot(4,2,5:6)
plot_confidence_intervals(x_axis/1000,ALL_M_norm(ALL_Hz_G == 10,:),[],[1,0,0])
plot_confidence_intervals(x_axis/1000,ALL_M_norm(ALL_Hz_G == 20,:),[],[0,1,0])
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
    plot_confidence_intervals(x_axis/1000,ALL_M_norm(ALL_LV_G == LVs(ii),:),[],clrs(ii,:))
end
plot_vert_line_at_zero
plot_vert_line_at_zero(-10);

legend_color_text({'0' '.3' '1' '1.2'},clrs(1:4,:))
title('LV')
xlabel('s')
ylabel('sd from base')

%% Plot the CV conditions

uLV = unique(ALLCV.LV);
uHz = unique(ALLCV.Hz);
for iLV = 1:length(uLV)
    for iHz = 1:length(uHz)
        Combo_IX{iLV,iHz} = ALLCV.Hz == uHz(iHz) & ALLCV.LV == uLV(iLV);
    end
end
figure
lab = [];
a = [];
clf
for iHz = 1:length(uHz)
    subplot(1,length(uHz),iHz)
    M = [];
    cnt = 1;
    for iLV = 1:length(uLV)
        RIX = Combo_IX{iLV,iHz}';
        if sum(RIX) > 1
            M(iLV,:) = mean(ALLCV.M(RIX,:));
            lab{iHz}{cnt} = sprintf('%1.2fLV',uLV(iLV));
            cnt = cnt + 1;
        end
    end
    plot(ALLCV.x_sec,M')
    xlabel('sec')
    title(sprintf('%s %1.2fHz', set_name,uHz(iHz)))
    axis tight
    a(iHz,:) = axis;

end
% equalize_axes;
aa = [a(1,1:2) min(a(:,3)) max(a(:,4))]
for iHz = 1:length(uHz)
    subplot(1,length(uHz),iHz)
    axis(aa)
    legend(lab{iHz}); legend boxoff;
end

%%
sm_bin = 10;
[pc,sc ] = pca(Z_scores(Q));
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
S.RT = movmean(S.RT,sm_bin,'omitnan');
S.Prop_Active = movmean(S.Prop_Active,sm_bin,'omitnan');
S.kurtosis = movmean(S.kurtosis,sm_bin,'omitnan');

win_bin = (20*1000)/binsize_ms;
warning off
[SD]= Sliding_dynamics_simple(Z_scores(Q), win_bin, round(win_bin/4));
warning on
SD_bin = [SD.win_bins(:,1) + round(win_bin/4)/2]';
bin_center_ms = round(win_bin/4)/2 * bin_ms;

SD_uS = linspace(Q_uS(1) + bin_center_ms*1000, Q_uS(end) - bin_center_ms*1000, length(SD.nEffDim));
nEffDim = interp1(SD_uS,SD.nEffDim,Q_uS,'linear');


figure
ax = [];
ax(1) = subplot(4,1,1);
imagesc(Q_uS/3600e6,[],Qz')
clim([-1 3])
xlabel('Hours'); ylabel('Neuron')
hold on
plot_markers_simple(TI.start_times_sec(TI.Hz_group == 10)/3600,[],1,[1,0,0])
plot_markers_simple(TI.start_times_sec( TI.Hz_group == 20)/3600,[],1,[0,1,0])
plot_markers_simple(TI.start_times_sec(TI.Hz_group == 60)/3600,[],1,[1,1,1])

plot_markers_simple(TI.end_times_sec(TI.Hz_group == 10)/3600,[],1,[1,0,0])
plot_markers_simple(TI.end_times_sec( TI.Hz_group == 20)/3600,[],1,[0,1,0])
plot_markers_simple(TI.end_times_sec(TI.Hz_group == 60)/3600,[],1,[1,1,1])

text(TI.end_times_sec/3600, Cols(Qz)*ones(1,length(TI.end_times_sec)), num2str(TI.LV_group))
text(TI.end_times_sec/3600, zeros(1,length(TI.end_times_sec)), num2str(TI.Hz_group))

yyaxis right
plot(cv_data(:,1)/3600,cv_data(:,2),'w.')
ylabel('[DA]')
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

figure
hold on
plot(cv_data(:,1)/3600,cv_data(:,2),'k.')
plot(Q_uS/3600e6,S.kurtosis)
yyaxis right
% plot(Q_uS/3600e6,S.RT)
plot(Q_uS/3600e6,pop_rate)

%%
figure
subplot(2,2,1:2)
plot(pc(:,1));hold on;plot(pc(:,2))
title('PC 1 an 2')
subplot(2,2,3)
plot(sc(:,1),sc(:,2),'.')
subplot(2,2,4)
plot(sc(:,2),sc(:,3),'.')
nEffDim(isnan(nEffDim)) = 0;
M= [Q_uS(:)/1e6 sc(:,1:2) pop_rate S.CV nEffDim];
sec_after = 280;
[TBL, INFO] = DANA_Condition_Aligned_Contin_Table_From_TI(M,{'PC1' 'PC2' 'pop_rate' 'CV' 'nEffDim'}, TI, rat_num, ALLCV, 'sec_after', sec_after);
%% This is great.
 data_type = 'PC1';
 data_type = 'PC2';
  % data_type = 'CV';
   % data_type = 'pop_rate';
  % data_type = 'nEffDim';
IX = TBL.Hz == 20 & TBL.LV ==1 & TBL.data_type == data_type; tstr = 'TBL.Hz == 20 & TBL.LV ==1';
IX = TBL.Hz > 10 & TBL.LV >=1 & TBL.data_type == data_type; tstr = 'IX = TBL.Hz > 10 & TBL.LV >=1';
 % IX = TBL.Hz > 0 & TBL.LV ==0 & TBL.data_type == data_type; tstr = 'IX = TBL.Hz > 0 & TBL.LV ==0';
% IX =  TBL.data_type == data_type & TBL.TrialNum == 1;
% IX =  TBL.data_type == data_type;

% TBL.data(IX,:)
ix = find(IX);
V = TBL.data(IX,:);
NIX = INFO.x_axis_s < -11;
Vn = V - mean(V(:,NIX),2);
Fx = INFO.FSCV_x_sec-15;

figure
subplot(4,1,1)
plot_raster(TBL.Stim_sequence_sec(ix)')
title(tstr)
pubify_figure_axis
box off 
axis ij

subplot(4,1,2)
imagesc(Fx, [], TBL.FSCV(IX,:))
set(gca,'XLim',INFO.x_axis_s([1 end]))
plot_vert_line_at_zero;plot_vert_line_at_zero(-10)
pubify_figure_axis
box off

yyaxis right
plot_confidence_intervals(Fx, TBL.FSCV(IX,:),[],[1 .8 .8])
set(gca,'XLim',INFO.x_axis_s([1 end]))
pubify_figure_axis
box off


subplot(4,1,3)
imagesc(INFO.x_axis_s ,[], V)
plot_vert_line_at_zero; plot_vert_line_at_zero(-10)
title(data_type)
yyaxis right
plot_confidence_intervals(INFO.x_axis_s, V)
pubify_figure_axis
ylabel(data_type)


subplot(4,1,4)
imagesc(INFO.x_axis_s ,[], Vn)
plot_vert_line_at_zero; plot_vert_line_at_zero(-10)
ylabel('Trial')

yyaxis right
plot_confidence_intervals(INFO.x_axis_s, Vn)
xlabel('sec')
pubify_figure_axis
set(gcf,'Position',[502    56   730   822])

%%

% Sparseness_measures(Q);

%% Correlations
window_size_bins = 20;
[self_similarity, out_times ]= Sliding_r(Q, Q_uS/3600e6, window_size_bins);
figure
plot(mean(out_times,2), self_similarity)


% % Local variance
% for iN = 1:length(T_uS)
%     LV(iN,:) = LocalVariance(diff(T_uS{iN},))
% end

%% END



if 0



    data_dir = fullfile(Git_dir,'\DANA\Data\Acute\20200527\Rat_B_unknown\52722_DANA_7000uM_bank0_g0');
    CV_long_csv_file = 'CV_data.csv';
    [NP, SP, CV] = Q2_does_DA_affect_ensemble_act_load(data_dir, CV_long_csv_file);
    binsize_ms = 20; % note - for analyses, we should really ignore bins that overlap the scans since those bins will always be zeros.
    sliding_window_ms = 4000; % for the dyanmics analysis and smoothing.
    PLOT_IT = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate the Q spike x time matrix.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_uS = []; All_T = [];
    for ii = 1:length(SP)
        T_uS{ii} = SP(ii).t_uS;
        All_T = [All_T; SP(ii).t_uS(:)];
    end
    All_T = unique(All_T);
    [a1,x1] = AutoCorr(All_T/10,1,800);
    [a2,x2] = AutoCorr(All_T/10,10,100);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;plot(x1,a1,x2,a2);xlabel('ms')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % There is a very strange 3 ms gap between samples - look at AutoCorr
    % No idea why. Need to figure this out. 10ms binning gets rid of the
    % effect, but why?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    depth_uM = [SP.depth_uM];
    [~,Qsix] = sort([SP.depth_uM]);
    depth_sorted_uM = depth_uM(Qsix);

    [Q_orig, edges_uS] = histcounts_cowen(T_uS, 'binsize', binsize_ms*1000);
    Q_uS = edges_uS(1:end-1) + (edges_uS(2) - edges_uS(1))/2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find bins that include scans (or stim - TODO:) and make them nans
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    edges_sec = edges_uS/1e6;
    BIX = false(size(edges_sec));
    E = [NP.scan_pulse_times_sec NP.scan_pulse_times_down_sec + 0.002];
    for ii = 1:length(NP.scan_pulse_times_sec)
        BIX(edges_sec >= E(ii,1) & edges_sec <= E(ii,2) ) = true;
    end

    QBIX = BIX(1:end-1);

    Q = Q_orig;
    Q(QBIX,:) = nan; % for the bins that overlap with the scans, make them nans
    Q_mm = movmean(Q,3,'omitnan' );

    figure
    a1 = subplot(3,1,1);
    imagesc(Q_orig'); colorbar
    caxis([0 3])
    a2 = subplot(3,1,2);
    imagesc(Q');
    caxis([0 3]); colorbar
    a3 = subplot(3,1,3);
    imagesc(Q_mm');colorbar
    caxis([0 3])

    linkaxes([a1 a2 a3])


    figure
    plot(mean(Q_orig,2))
    hold on
    plot(mean(Q,2))
    plot(mean(Q_mm,2))


    [xc1,~] = xcorr(mean(Q_orig,2),100);
    [xc2,~] = xcorr(mean(Q,2),100);
    [xc3,lags] = xcorr(mean(Q_mm,2),100);

    figure
    plot(lags*binsize_ms,xc1)
    hold on
    plot(lags*binsize_ms,xc2)
    plot(lags*binsize_ms,xc3)
    xlabel('ms lag')
    legend('orig', 'with nan', 'movmean')

    % Determine higher-order states.
    [~,Qpc,Qlat] = pca(Q_mm); % I am a little suspicious as like 95+% of data is concetrated in the first component.
    % [~,Qpc_tmp,Qlat] = pca(Q_mm);
    % Qpc = nan(size(Q));
    % Qpc(~QBIX,:)= Qpc_tmp;
    window_size_bins = sliding_window_ms/binsize_ms;
    SD = Sliding_dynamics_simple(Q(~QBIX,:), window_size_bins);
    SD.win_bin_centers_ix = round(mean(SD.win_bins(),2));
    good_sec = Q_uS(~QBIX)/1e6;
    SD.win_bin_centers_sec = good_sec(SD.win_bin_centers_ix);
    % IFR
    [IFR,IFRsmth,bin_edges_s] = Instantaneous_Firing_Rate(NP.stim_times_sec,.1,5);
    IFR_sec = bin_edges_s(1:end-1) + (bin_edges_s(2)- bin_edges_s(1))/2;

    if PLOT_IT
        figure
        a1 = subplot(3,1,1);
        imagesc(Q_uS/60e6,[],Q(:,Qsix)')
        caxis([0 2])

        a2 = subplot(3,1,2);
        plot(SD.win_bin_centers_sec/60,SD.nEffDim)
        yyaxis right
        plot(SD.win_bin_centers_sec/60,SD.R_matrix_mn_r)
        axis tight
        title('nEffDim, Rmn r')


        a3 = subplot(3,1,3);
        plot(SD.win_bin_centers_sec/60,SD.CV)
        ylabel('CV')
        yyaxis right
        plot(SD.win_bin_centers_sec/60,SD.prop_active)
        ylabel('prop active')
        axis tight

        %     linkaxes([a1 a2 a3]) % does not work with yyaxis

        % something seems a little screwy with IFR - does not map onto what I see sometimes.
        figure
        plot(IFR_sec,IFRsmth,IFR_sec,IFR)
        hold on
        plot(NP.stim_times_sec,zeros(size(NP.stim_times_sec)),'r+')


        figure
        plot(CV.NP_times_sec,CV.Concentration)

        % Map this onto the Q data using interp 1.

        figure
        ax1 = subplot(3,1,1:2);
        imagesc(Q_uS/60e6,[],Q_mm(:,Qsix)');
        caxis([0 5])
        ylabel('neuron ID (sort by depth)')
        yyaxis right
        colororder({'k','b'})
        plot(CV.NP_times_sec/60,CV.Concentration,'w','LineWidth',2)
        xlabel('min')
        ylabel('[DA]')
        pubify_figure_axis
        ax2 = subplot(3,1,3);
        plot(IFR_sec/60, IFRsmth,'-','LineWidth',2,'Color',[.8 .1 .1])
        hold on
        plot(NP.stim_times_sec/60,zeros(size(NP.stim_times_sec)),'r+')
        xlabel('min')
        ylabel('IFR')
        yyaxis right
        colororder({'k','b'})

        plot(Q_uS/60e6,mean(Q_mm,2),'k-','LineWidth',1)
        %      hold on
        %      plot(Q_uS/60e6,movmedian(Qpc(:,1),10),'b-')
        ylabel('mean activity')
        axis tight
        pubify_figure_axis
        linkaxes([ax1,ax2],'x')
    end
    %%
    % Now to things by condition.
    uLV = unique(NP.LV(NP.Hz == 20));
    FSCV_DATA = [CV.NP_times_sec,CV.Concentration];
    GIX = ~isnan(sum(FSCV_DATA,2));
    FSCV_DATA = FSCV_DATA(GIX,:);
    GIX = ~isnan(sum(Qpc,2));
    PC_DATA = [Q_uS(GIX)/1e6 movmedian(Qpc(GIX,1),50)];
    Q_DATA = [Q_uS/1e6 mean(Q_mm,2)];
    IF_DATA = [IFR_sec IFR IFRsmth];
    sFreq_IFR = 1/median(diff(IFR_sec));
    % smooth the PC a bit.
    sFreq_FSCV = 5;
    sFreq_Q = 1/median(diff(Q_uS(GIX)/1e6));
    sec_before = 5; sec_after = 20;

    for ii = 1:length(uLV)
        IX = NP.LV ==uLV(ii) & NP.Hz == 20;
        [FSCV, ~, x_sec_FSCV] = PETH_EEG_simple(FSCV_DATA, NP.Stim_Start_sec(IX), ...
            sec_before*sFreq_FSCV, sec_after*sFreq_FSCV, sFreq_FSCV,false);
        [PC, ~, x_sec_PC] = PETH_EEG_simple(PC_DATA, NP.Stim_Start_sec(IX), ...
            sec_before*sFreq_Q, sec_after*sFreq_Q, sFreq_Q,false);
        [MNQ, ~, x_sec_Q] = PETH_EEG_simple(Q_DATA, NP.Stim_Start_sec(IX), ...
            sec_before*sFreq_Q, sec_after*sFreq_Q, sFreq_Q,false);
        [IF, ~, x_sec_IF] = PETH_EEG_simple(IF_DATA(:,1:2), NP.Stim_Start_sec(IX), ...
            sec_before*sFreq_IFR, sec_after*sFreq_IFR, sFreq_IFR,false);
        [IFs, ~, x_sec_IF] = PETH_EEG_simple(IF_DATA(:,[1 3]), NP.Stim_Start_sec(IX), ...
            sec_before*sFreq_IFR, sec_after*sFreq_IFR, sFreq_IFR,false);

        figure
        subplot(5,1,1)
        imagesc(x_sec_FSCV,[],FSCV)
        ylabel('trial')
        pubify_figure_axis
        title(sprintf('LV = %1.2f, 20Hz stim, FSCV',uLV(ii)))

        subplot(5,1,2)
        plot_confidence_intervals(x_sec_FSCV,FSCV);
        ylabel('[DA]')
        pubify_figure_axis

        subplot(5,1,3)
        imagesc(x_sec_PC,[],PC)
        pubify_figure_axis
        title('PC1 of neural ensemble')

        subplot(5,1,4)
        plot_confidence_intervals(x_sec_Q,MNQ,[],[.1 .2 .9]);
        ylabel('mean activity')
        yyaxis right
        plot_confidence_intervals(x_sec_PC,PC);
        ylabel('PC1')
        xlabel('sec')
        pubify_figure_axis

        subplot(5,1,5)
        plot_confidence_intervals(x_sec_IF,IF);
        hold on
        plot(x_sec_IF, mean(IFs),'r','LineWidth',2)

        ylabel('IFR')
        xlabel('sec')
        pubify_figure_axis

        set(gcf,'Position',[489 49 560 821.6])
    end


    nEffDim_DATA = [SD.win_bin_centers_sec SD.nEffDim];
    sFreq_nEffDim = 1/median(diff(SD.win_bin_centers_sec));

    for ii = 1:length(uLV)
        IX = NP.LV ==uLV(ii) & NP.Hz == 20;
        [FSCV, ~, x_sec_FSCV] = PETH_EEG_simple(FSCV_DATA, NP.Stim_Start_sec(IX), ...
            sec_before*sFreq_FSCV, sec_after*sFreq_FSCV, sFreq_FSCV,false);
        [nEffDim, ~, x_sec_nEffDim] = PETH_EEG_simple(nEffDim_DATA, NP.Stim_Start_sec(IX), ...
            sec_before*sFreq_nEffDim, sec_after*sFreq_nEffDim, sFreq_nEffDim,false);


        figure
        subplot(5,1,1)
        imagesc(x_sec_FSCV,[],FSCV)
        ylabel('trial')
        pubify_figure_axis
        title(sprintf('LV = %1.2f, 20Hz stim, FSCV',uLV(ii)))

        subplot(5,1,2)
        plot_confidence_intervals(x_sec_FSCV,FSCV);
        ylabel('[DA]')
        pubify_figure_axis

        subplot(5,1,3)
        imagesc(x_sec_nEffDim,[],nEffDim)
        pubify_figure_axis
        title('PC1 of neural ensemble')

        subplot(5,1,4)
        plot_confidence_intervals(x_sec_nEffDim,nEffDim);
        ylabel('nEffDim')
        xlabel('sec')
        pubify_figure_axis

        subplot(5,1,5)
        plot_confidence_intervals(x_sec_IF,IF);
        hold on
        plot(x_sec_IF, mean(IFs),'r','LineWidth',2)

        ylabel('IFR')
        xlabel('sec')
        pubify_figure_axis

        set(gcf,'Position',[489 49 560 821.6])
    end

    % INFO = PETH_plot_and_analyze_population(Q_uS/1e6,Qz')
end