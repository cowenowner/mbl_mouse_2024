%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HI!!
% Here is a demonstration that will teach you a little about analyzing nerual
% ensemble data.
%
% EXPERIMENT BACKGROUND:
% You can check out the poster in this folder, but here is the short
% version: Electrical brain stimulation of neural circuits connected to the
% reinforcement system (e.g., VTA, MFB, NAc) has been used for over 70
% years to motivate animals to press levers and brain stimulation of the
% STN and GPi has been used for decades to treat Parkinson's disease.
% Almost all of the standard procedures for stimulating the brain involve
% delivering predictable pulse trains with identical inter-pulse intervals:
% |   |   |   |   |   |
%
% The issue is that neural circuits and individual neurons are adaptable
% and can habituate to predictable stimuli. Why wouldn't this be the case
% for brain stimulation? We asked that question but more specifically, we
% asked wheter variable inter-pulse intervals with statistics that resemble
% 'bursts' would be the most effective at evoking dopamine release and
% neural ensemble activity, even if the mean stimulation frequeny is held
% constant. For example...
%
% |      ||  |  |     |
%
% How to quantify inter-pulse variability? There are many ways. For
% example, you could just calculate the standard deviation of the
% inter-pulse intervals. Other measures include the coefficient of
% variation or the Fano factor. I chose 'local variance' (Shinomoto et al.)
% which measures inter-pulse variabilty but controls for slower changes in
% firing rate. LV of 0 indicates no inter-pulse variance, LV of 1 = Poisson
% inter-pulse variance, LV > 1 indicates 'bursty' inter pulse variance
% (many shorter intervals but some LONG inter pulse intervals.
%
% BRAIN REGIONS: In this experiment, rats are anesthetized and we have a
% bipolar stimulation electrode implanted in the medial forebrain bundle
% (MFB) and a neuropixels probe implanted in the dorsal and ventral
% striatum (lower channel numbers = deeper electrodes). The ventral
% striatum is also called the nuclues accumbens (NAc). In some animals, we
% were also able to simultaneously measure dopamine using fast-scan-cyclic
% voltammetry (FSCV). We'll save that for another demo however.
%
% GOALS AND TASKS:
%
%        1) Get the code to run once by calling 
%           DEMO1_does_MFB_stim_affect_ensemble_act on the command line. 
%           This will make sure it all runs before we dig into the demo.
%        2) Create a folder in mbl_mouse_2023\student_code with your name.
%        3) In that folder, create a powerpoint file called
%           'Analysis_notebook_for_DEMO1.ppt'. This file will be useful for
%           keeping your insights and notes regarding your analyses.
%        4) Copy this demo to that directory, but put your name at the end.
%           DEMO2_does_MFB_stim_affect_ensemble_act_cowen.m
%        5) Now step through the code, line by line, and try to get into
%           the mind of the crazy person who wrote the code. Since this is now
%           your copy of the code, feel free to mess around with stuff. Try
%           new things. I would just copy a few lines at a time and run them
%           in the CommandWindow.
%        6) As you go through the code, answer each question posed in the
%           code (indicated with a ?), but PLEASE ask your friends,
%           neighbors, assistants, sketchy randos, or me if you have any 
%           questions. I am happy to explain. I also may randomly quiz 
%           you!! When you least expect it, expect it! Also, don't forget 
%           that you can type doc <matlabfunction> to get the documentation 
%           for that function. There is always ChatGPT (which is
%           suprisingly good at matlab)
%
% An important thing to know: The method we used to collect dopamine (FSCV)
% requires sending a very large pulse (1.5Volts) into the brain every 200
% ms (pulse width is 8ms). This creates a large artifact that prevents
% neural recording at the time. Similarly, each MFB stimulation creates a
% similar artifact.
%
% DEPENDENCIES:  (on GitHub)
% This code requires my lab's communal code repository: CowenLib.
% It also requires the code specific to the DANA project: DANA/src/matlab
% (add paths and all subfolders).
%
% CONVENTIONS:
% General purpose functions that my lab has made or borrowed from others
% are in CowenLib. Functions that are specific to a particular project or
% experiment typically have a prefix in all caps such as
% DANA_trial_info_from_stim_times.m. This makes it easy to separate general
% utility funcitons from experiment specific functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The name of the data file to load. There are a few data files in this
% directory. Feel free to try them all and look at the differences.
clearvars
filename = 'DANA_rat445_01.mat'; % Assumes that this file is in the directory with the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters to play with...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ? What is a vector, matrix, cell array, structure, boolean?
big_binsize_ms = 500; % Binsize for spike x time matrices.
small_binsize_ms = 20; % Binsize for spike x time matrices.
% ? What do you think would be a good choice of bin size.
PLOT_IT = true; % If you don't want to plot and just want to process 
selfsim_win_size_bins = 20; % for an esoteric analysis at the end.
plot_each_peth = false; % boolean for controlling plots. Sometimes too many plots are just too many.
time_before_sec = 25; % time before and after the peri-event-time histogram (PETH)
time_after_sec = 60*2; % see above.
peth_bin_ms = 550;
clrs = lines(200); % this stores some nice colors to use for plotting the spike rasters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data.
load(filename);
% Put timestamps into a nice format. SP: a structure array (very similar to
% a database) that stores each cell's timestamps, waveformshape, auto-corr,
% and a host of other things.

% EVT: a structure that has the times of stimulation pulses and the FSCV
% scan pulses. ?What is FSCV?
%%
T_uS = {SP.t_uS}; % This converts the timestamp list for each cell into a more easy to use cell array.
WV = [SP.WaveformBest]'; % Waveform shapes for each neuron. This could indicate cell type.

% Now that the data is loaded, make sure that it looks pretty. This uses my
% plot_raster() function. Feel free to use it if you like it.
figure; plot_raster(T_uS,[],'time_divisor',60e6);xlabel('minutes')
%%
%%%%%%%%%%%%%%%%%
% Consider each 10-second stim pulse train a 'trial' as it has a unique set
% of features (stim frequency, local variance). Conditions are interleaved.
% Randomly interleaving trial conditions is almost always the best idea for
% many reasons such as helping deal with slow 'non-stationarities' in the
% data.
% ? What is an example of a non-stationarity in neural data?
% Determine trial information (e.g., stim frequency, local variance) from the stim times.
% Information includes start/end time, stim freq, etc...
% ? What is a structure data type in matlab?
%%%%%%%%%%%%%%%%%
% TI = DANA_trial_info_from_stim_times(stim_times_sec); % Feel free to look at this function - only if you feel brave though. It's dense.
% The TI structure has TI.LV (the un-normalized local variance) and TI.LVR
% (an improved measure of local variance from a more recent paper).
%%
% Create a neuron x time matrix. Let's call it 'Q'
[Q, edges_uS] = histcounts_cowen(T_uS, 'binsize', big_binsize_ms*1000);
Q_uS = edges_uS(1:end-1) + (edges_uS(2) - edges_uS(1))/2;

size(Q)

figure;
subplot(3,1,1);
imagesc(Q_uS/60e6,[],Q')
ylabel('Neuron');xlabel('time (min)')
colorbar
title(sprintf('Original data binned at %d ms', big_binsize_ms))
% ?What did the ' do?
% ?Why did I add /60e6? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's normalize each neuron's activity.
% Normalization or standardization can be done in many ways. Consider
% carefully what the implications of normalization might be for
% interpreting your data. For example, log normalization reduces the impact
% of very high firing rates. If you think that very high firing rates are
% very informative, then this might not be a good idea. z scores assume
% that the absolute mean rate of a cell is not important as the mean is
% subtraced.
% Z score: subtract mean and divide by standard deviation so that data is in units
% of standard deviation from the mean
%
% Log: take the logarithm of the data.
Qz = Z_scores(Q);
subplot(3,1,2);
imagesc(Q_uS/60e6,[],Qz')
colorbar
ylabel('Neuron');xlabel('time (min)')
title('z score standardization')

Qlog = log10(Q+1); 
% ?Why do I need to add 1?
subplot(3,1,3);
imagesc(Q_uS/60e6,[],Qlog')
colorbar
ylabel('Neuron');xlabel('time (min)')
title('log norm')

%% Let's get a visualization of network-level interactions in the form of 
% correlations between cell pairs.
R = corr(Q); 
% ?What are the units of R?
figure
subplot(1,2,1)
imagesc(R)
xlabel('Neuron');ylabel('Neuron')
clim([-.1 .5]) % allows you to control the color scale
colorbar
% Now let's look at how this changes when we choose a smaller time bin.
[Qsmallbin, small_edges_uS] = histcounts_cowen(T_uS, 'binsize', small_binsize_ms*1000);
Q_small_uS = small_edges_uS(1:end-1) + (small_edges_uS(2) - small_edges_uS(1))/2;

Rsmall = corr(Qsmallbin);
subplot(1,2,2)
imagesc(Rsmall)
xlabel('Neuron');ylabel('Neuron')
clim([-.1 .5])
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's get back to the experiment and plot the data by condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This plot summarizes the entire experiment by marking each stimulatoin
% episode and the rate and LV for that stim.

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
%%

% Create table of aligned neural responses. This function does the hard work of 
% creating peri-event time histograms. It's a rather complicated function
% so you might just have to trust it.
% [TBL, INFO] = DANA_Condition_Aligned_Spike_Table_From_TI(SP,TI,rat_num,'time_after_sec',time_after_sec);
% In addition to matrices, cell arrays, and structures, matlab also has
% 'tables' which are similar to excel spreadsheets or R data frames.
TBL
% look at the output of TBL - that will help you understand what it stores.
% The following plot asks: Does stimulation alone or frequency of
% stimulation have an effect on the firing rate during stimulation or after
% stimulation?
% Stimulation off is at 0 seconds.
tmp_clrs = [1 .2 .2;.2 .2 1]; % Colors RGB
figure
uHz = unique(TBL.Hz);
for ii = 1:length(uHz)
    IX = TBL.Hz == uHz(ii);
    plot_confidence_intervals(INFO.PETH_x_axis_ms/1000,TBL.mean_PETHnorm(IX,:),[],tmp_clrs(ii,:))
end
title('Hz'); ylabel('Norm rate'); xlabel('sec')
pubify_figure_axis
legend_color_text({'10 Hz' '20 Hz'},tmp_clrs(1:2,:))
% OK, a bit interesting. Somehow 10Hz results in a bigger rebound after
% stim offset.

% Now let's ask the same question but for local variance of the
% stimulation...

figure
uLVs = unique(TBL.LV);
for ii = 1:length(uLVs)
    IX = TBL.LV == uLVs(ii);
    plot_confidence_intervals(INFO.PETH_x_axis_ms/1000,TBL.mean_PETHnorm(IX,:),[],clrs(ii,:))
end
title('LV'); ylabel('Norm rate'); xlabel('sec')
pubify_figure_axis
legend_color_text({'0 LV' '.3' '1' '1.2'},clrs(1:4,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Network measures part 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The above was a single-cell analysis in the sense that we computed a
% measure for each neuron (a peri-event response) and then just averaged
% across neurons and by group (e.g., 10Hz and 20Hz). Since one theme for
% this mouse cycle will be network-level analyses and brain-state
% identification, let's get a taste for network level measures. We already
% started this with the R matrix - a very useful visualization of
% interactions between cells, but there are many other measure of
% interactions.

% Perhaps the most well known measure of network activity (a network of 2
% cells so a dinky network) is the cross-correlation.
% cell_1 = 1; cell_2 = 2;
cell_1 = 3; cell_2 = 4;
xc_x_max_ms = 200;
xc_bin_ms = 5;
n_lags = round(xc_x_max_ms/xc_bin_ms);
[cc,x] = Cross_corr(T_uS{cell_1}/1e3,T_uS{cell_2}/1e3,xc_bin_ms,n_lags);

figure
plot(x,cc)
plot_vert_line_at_zero
xlabel('ms')
title('Cross Correlogram')
ylabel('n coincident events')

% Some fun explorations:
% 1) now try some different cell combinations
% 2) now make a few auto-correlograms
% 3) For the brave, write a for loop (or a for loop within a for loop) that
% does every possible combination of cross correlograms.

%% PCA
% The cross-corr is the workhorse but its weakness is that it can only look at 2 cells at a time. 
% More modern approaches use different approaches for either 1) reducing
% the dimensionality of the Q matrix above (from n cells to something much
% lower dimension) or 2) cluster the population vectors (slices of the Q
% matrix) into categorcical/discrete states. The first is often done by
% pca,ica,auto-encoders. The second is often done by kmeans, hierarchical
% clustering, or more fancy things that I don't understand.

% Let's start with PCA...
% Guaranteed the following paragraph will be confusing - so please google
% PCA to get a better description.
%
% PCA identifies ways to 'rotate' your high dimensional data into lower
% dimensional spaces. Another way to think of PCA is that it finds new ways
% to compress your high dimensional (n neuron dimensions) into a lower
% dimensional representation - like compression of music to an mp3. 
% For a time x neuron matrix, PCA finds n-neurons ways to compress your
% data. It ranks those n-ways by how well they compress the data. You can
% see the degree of compression in the 'lat' output in the table below.
% Each column of pc indicates a template or 'eigenvector' that is applied
% via dot product to each time slice of your neuron x time matrix to create
% the new dimension. The sc matrix is a new matrix with each column
% representing the lower dimensional version of your neuron x time matrix.
[pc,sc,lat] = pca(Qsmallbin);
figure
subplot(1,3,1);imagesc(pc);xlabel('pc num')
subplot(1,3,2);imagesc(sc);ylabel('time')
subplot(1,3,3);bar(lat); ylabel('variance explained')
% let's plot the 'best' lower dimensional representaion of the neuron x
% time matrix...
figure
plot(Q_small_uS/60e6,sc(:,1))
xlabel('min')
% now let's plot the second best.
hold on
plot(Q_small_uS/60e6,sc(:,2))
legend('PC1','PC2')
%% Cluster plots
% You can also visually check to see if there are 'clusters' of states in
% the data by making a scatter plot of the first few dimensions. 
figure
subplot(2,2,1)
plot(sc(:,1),sc(:,2),'.'); xlabel('PC1');ylabel('PC2')
subplot(2,2,2)
plot(sc(:,2),sc(:,3),'.'); xlabel('PC2');ylabel('PC3')
subplot(2,2,3)
plot(sc(:,3),sc(:,4),'.'); xlabel('PC3');ylabel('PC4')
subplot(2,2,4)
plot(sc(:,1),sc(:,3),'.'); xlabel('PC1');ylabel('PC3')
% These look crazy BTW.

% for kicks, try the above with the big binned Q matrix (500 ms bin)

%% Now suppose we would like to identify 'discete' or categorical states or 'clusters' in this high-dimensional data. How could we do this? There are MANY ways
% let's start with the simplest or at least the most common: kmeans.
% Since I will murder the explanation, please google kmeans or ask ChatGPT.
n_clusters = 5;
C = kmeans(Qsmallbin,n_clusters);
% C indicates the cluster identity for each time slice (population vector)
% in the Q matrix. It is a number from 1:n_clusters
figure
histogram(C)
% ?are there more of a certain type of cluster?
% Let's dig deeper into the clusters. What do they look like?
mean_pop_vec = [];
for iK = 1:n_clusters
    IX = C == iK;
    mean_pop_vec(iK,:) = mean(Qsmallbin(IX,:));
end
% FYI a far more 'compact' way to do the above would be...
mean_pop_vec2 = grpstats(Qsmallbin,C,'mean');
% but it's good to learn for loops as they can generalize to many other
% situations.
figure
subplot(2,1,1)
imagesc(mean_pop_vec')
xlabel('Cluster number')
subplot(2,1,2)
bar(mean(mean_pop_vec,2))
axis tight
ylabel('mean rate')

% Let's put our clusters back in PC space...
figure
gscatter(sc(:,1),sc(:,2),C); xlabel('PC1');ylabel('PC2')

%
% I think that the first state is the most common according to the histogram. 
% Is there anything special about this state?
% No try again with the big bins. Do you get a more equal distribution of
% states?
%
%% Now let's look at these states over time to see if there is a pattern
% Feel free to zoom in and explore! Is there a pattern to the states?
% Heck if I know - how would you test this?
ncells = size(Qsmallbin,2);
figure
a = [];
a(1) = subplot(4,1,1:3);
% ax(1) = subplot(4,1,1);
imagesc(Q_small_uS/3600e6,[],Qsmallbin')
clim([-1 3])
xlabel('Hours'); ylabel('Neuron')
hold on
plot(stim_times_sec/3600,ones(size(stim_times_sec)),'w^')
plot_markers_simple(TI.end_times_sec(TI.Hz_group == 10)/3600,[],1,[1,0,0])
plot_markers_simple(TI.end_times_sec( TI.Hz_group == 20)/3600,[],1,[0,1,0])
plot_markers_simple(TI.end_times_sec(TI.Hz_group == 60)/3600,[],1,[1,1,1])

text(TI.end_times_sec/3600, ones(1,length(TI.end_times_sec)), num2str(TI.LV_group))
text(TI.end_times_sec/3600, zeros(1,length(TI.end_times_sec)), num2str(TI.Hz_group))
sgtitle(sprintf('%s Population activity (R = 10Hz, G = 20Hz)',set_name))

a(2) = subplot(4,1,4);
hold on
% plot(Q_small_uS([1 end])/3600e6, [ncells+1 ncells+1],'w','LineWidth',10)
for iK = 1:n_clusters
    IX = C == iK;
    times = Q_small_uS(IX)/3600e6;
    plot(times, zeros(size(times))+iK,'.','MarkerSize',23,'Color',clrs(iK,:))
    plot(times, zeros(size(times))+iK,'.','MarkerSize',6,'Color',[1 1 1])

end
axis tight
linkaxes(a,'x')
