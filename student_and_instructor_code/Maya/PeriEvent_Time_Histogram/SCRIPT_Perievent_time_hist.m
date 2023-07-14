%% Peri-Event Time Histograms 
% This script is used to create histograms locked on a particular
% event/stimulus presentation. Code is from MVMD
clearvars; clear path; clear classes
%% Load the code
% addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\CowenLib'));
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\mvdm'));
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\Maya\PeriEvent_Time_Histogram'));
addpath(genpath('C:\data\02_7_7_23\preprocessed'));
%% Load the data
%clear;
datadir = ('C:\data\02_7_7_23\preprocessed'); %Create data directory

load(fullfile(datadir,'AllSpikes.mat')); % Load in experiment data
load(fullfile(datadir,'odor_events.mat')); %Load in odor presentation times
load(fullfile(datadir,'LFPs.mat')); %Load in LFP data HUUUGE


%% Defining Multi-Unit Activity (MUA)
S = ConvertSPtoS(SP); %Use Function to convert SP to S

cfg_MUA = []; %Create empty variable for multi-unit activity(MUA)
cfg_MUA.tvec = lfp_tsd.tvec'; % Synchronize timebase to compute MUA on
MUA = getMUA(cfg_MUA, S); %Make vector with MUA, S (?)

%% Comparing MUA waveform with z-transformed MUA waveform
%**How subplot 121 works:you would like to have one row and two columns worth
% of figures. The last number, p=1 means that you wish to place the plot in
% the left most column**

figure; %Make a figure
subplot(121) %Subplot the data
plot(MUA) %Plot MUA
title('raw MUA (spikes/s)') %Add title

subplot(122) %subplot more data!
MUAz = zscore_tsd(MUA); %Convert MUA into z-score values
plot(MUAz) %Plot the z score values
title('z-scored MUA');%Add title

%% Event-triggered MUA

for iC = 1:3 %For cue 1-3
    cue_evt(iC) = SelectTS([], evt, iC); % select events for each individual cue
end

cfg_peth = []; % Create vector containing parameters for PETH
cfg_peth.window = [-1 1]; %Specify window going from -1 to 1
cfg_peth.dt = 0.01; %Idunno
cfg_peth.mode = 'interp'; %Interpolate input data to smooth it out

figure; %Make a new figure
cols = 'rkb'; %Specify colors red, black and blue
for iC = 1:3 %For cue 1-3
    out = TSDpeth(cfg_peth, MUAz, cue_evt(iC)); %Make a vector "out" equal to the cfg_peth vector, MUAz, and cue times
    h(iC) = plot(out, 'Color', cols(iC), 'LineWidth', 2); %Plot vector "out", with color and linewidth
    hold on; %Hold on graph; loops for each cue 1-3
end

set(gca, 'FontSize', 18, 'TickDir', 'out'); box off; %Setting the axes properties
legend(h, {'Odor 1', 'Odor 2', 'Odor 3'}); legend boxoff; %Make a legend box

xlabel('Time from cue onset (s)') %Label x axis
ylabel('z-scored Multi-Unit Activity') %Label y axis

vline(0,':'); %Addd a line at x = 0

%% Single unit Peri-event histograms
%To plot 
for iCa = 10:18 % What does this line do?
    this_S = SelectTS([], S, iCa); % select only the current cell
    MUA = getMUA(cfg_MUA, this_S); % "MUA" for one cell is just that cell's firing rate
    MUAz = zscore_tsd(MUA); %Transform into zscore
    subplot(3,3,iCa) %Make a plot with subplots

    for iCue = 1:3 %For cue 1-3
        this_out = TSDpeth(cfg_peth, MUAz, cue_evt(iCue)); hold on; %make vector "this out"
        plot(this_out, 'Color', cols(iCue), 'LineWidth', 2); %plot "this out"
        title(sprintf('cell %d', iCa));
    end
    vline(0,':');
end