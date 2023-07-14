%% LFP/SWR Visualization
% This script is used to visualize the data in order to collect examples of
% SWRs in the dataset. The script can be adapted depending on where you
% want to survey SWRs (i.e. which shank, what subset of neurons). Code is 
% from MVMD
clearvars; clear path; clear classes
%% Load the code
%Data should be stored on Github repository
%difference between addpath and genpath?

addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));
addpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\mvdm');
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\Maya\SWR_Visualization'));

%% Load the Data
%Data should be stored on local folder, NOT GITHUB.
clear;
datadir = ('C:\data\02_7_7_23\preprocessed'); %Create data directory

load(fullfile(datadir,'AllSpikes.mat')); % Load in experiment data
load(fullfile(datadir,'odor_events.mat')); %Load in odor presentation times
load(fullfile(datadir,'LFPs.mat')); %Load in LFP data

%% Create event time vectors
%Create a vector containing stimulus presentation times
event_time1 = evt.t{1}; 
event_time2 = evt.t{2};
event_time3 = evt.t{3};

%% Check struct of variable lfp_tsd (from LFPs.mat)
%Contains lfp data with timestamps in seconds
%Output= type, units, tvec (time vector?), data(how many channels)
lfp_tsd

%% Plot coordinates to visualize neuropixel shanks
%Will show the shanks from which the dataset was recorded from, a total of
%384 sites divided into 4 (?)

plot(lfp_tsd.usr.xcoord, lfp_tsd.usr.ycoord, 'ok')

%% Subsetting data with keep.idx

keep_idx = lfp_tsd.usr.xcoord == 277; %Make index associated with the probe x-axis you want to visualize
lfp_tsd.data = lfp_tsd.data(keep_idx, :); %Make vector to synchronize the other datasets with this index
lfp_tsd.usr.xcoord = lfp_tsd.usr.xcoord(keep_idx); %Synchronizing x coordinates of probe
lfp_tsd.usr.ycoord = lfp_tsd.usr.ycoord(keep_idx); %Synchronizing y coordinates of probe
lfp_tsd  %Confirm by checking struct; should be 96 if only choosing one shank (?)


%% Sorting Data

plot(lfp_tsd.usr.ycoord) % Checking data channels are sorted in order
%If not a straight line, do the following sorting:

[~, sort_idx] = sort(lfp_tsd.usr.ycoord, 'ascend'); % Create function to sort values ascending
lfp_tsd.usr.xcoord = lfp_tsd.usr.xcoord(sort_idx); % Sort x coordinate data
lfp_tsd.usr.ycoord = lfp_tsd.usr.ycoord(sort_idx); % Sort y coordinate data
lfp_tsd.data = lfp_tsd.data(sort_idx, :); % Synchronize lfp_tsd.data with sorted values
plot(lfp_tsd.usr.ycoord) % Checking data channels are sorted in order

%% Thinning out visualized channels
% Using a subset of the LFP data so processing is faster

keep_idx = 1:3:96; % Creating index that includes only every third channel
lfp_tsd.data = lfp_tsd.data(keep_idx, :); % Synchronize lfp_tst.data with subset
lfp_tsd.usr.xcoord = lfp_tsd.usr.xcoord(keep_idx); %Synchronize x coordinate data with subset
lfp_tsd.usr.ycoord = lfp_tsd.usr.ycoord(keep_idx); %Synchronize y coordinate data with subset
lfp_tsd % Confirm by checking struct; should be 1/3rd of the total number LFPs

%% LFP Filtering

plot(lfp_tsd.data(17,:)) % Plotting a single LFP channel

chan_list = 1:32; %Creating variable with the channel list
%Loop to create new vector with filtered LFP

nCh = length(chan_list);
cfg_mr = [];
for iCh = 1:nCh; % For every channel in the channel list
    cfg_mr.lfp(iCh) = lfp_tsd; % Structure cfg_mr.lfp(iCh) equal to lfp_tsd variable
    cfg_mr.lfp(iCh).data = cfg_mr.lfp(iCh).data(chan_list(iCh), :); % Structure cfg_mr.lfp(iCh).data equal to all columns for each row 
    cfg_f = []; %Create an empty vector
    cfg_f.f = [1 475]; % filter in LFP range, this makes the plots easier to read
    cfg_mr.lfp(iCh) = FilterLFP(cfg_f,cfg_mr.lfp(iCh)); % Apply filter to data
    %cfg_mr.spkColor = 'w';
end

%If getting errors:
%CheckTSD(cfg_mr.lfp(iCh)) % Checks if variable is in the proper format (should return 1)
%% Spike timestamp: convert to ordered ts in seconds
S = ts; % Make variable myS equal to timestamp data

% Synchronizing time variable with cell clusters for LFPs
nCells = length(SP); %Define number of cells; found in SP struct length

% convert to s and align with LFPs
for iC = 1:nCells %For timestamps of each cell
    SP(iC).t = SP(iC).t_uS * 10^-6; %Multiply by 10^-6 to get seconds
end

% make timestamp struct for spikes using a loop
for iC = nCells:-1:1 %For each cell in list
    S.t{iC} = SP(iC).t; %Get time variable in S equal to time variable in SP
    S.label{iC} = SP(iC).cluster_id; %Make sure clusters in SP align with labels in S
end

%% Plot the raster
doc MultiRaster

 %Subset Data
 % Baseline = 

MultiRaster(cfg_mr, S);
set(gca, 'Color', [0 0 0]);
set(gca, 'XColor', [1 1 1]);
set(gca, 'YColor', [1 1 1]);
set(gcf, 'Color', [0 0 0]);
set(gca, 'YTick', [])
ylabel('')
hold on; 
arrayfun(@(x)xline (x, '-', '1', 'LabelOrientation', 'horizontal'), event_time1);%-- plot vertical lines for event1
arrayfun(@(x)xline (x, '--', '2', 'LabelOrientation', 'horizontal'), event_time2);%-- plot vertical lines for event2
arrayfun(@(x)xline (x, ':', '3', 'LabelOrientation', 'horizontal'), event_time3);%-- plot vertical lines for event3
%% Plotting around each event
this_event_id = 54;
this_ts = event_time2(this_event_id);
% restrict the tsd and S to some time around the event
this_lfp_tsd = restrict(lfp_tsd, iv(this_ts-2,this_ts+2));
this_S = restrict(S, iv(this_ts-2,this_ts+2));

cfg_mr2 = [];
for iCh = 1:nCh % For every channel in the channel list
    cfg_mr2.lfp(iCh) = this_lfp_tsd; % Structure cfg_mr.lfp(iCh) equal to lfp_tsd variable
    cfg_mr2.lfp(iCh).data = cfg_mr2.lfp(iCh).data(chan_list(iCh), :); % Structure cfg_mr.lfp(iCh).data equal to all columns for each row 
    cfg_f2 = []; %Create an empty vector
    cfg_f2.f = [1 475]; % filter in LFP range, this makes the plots easier to read
    cfg_mr2.lfp(iCh) = FilterLFP(cfg_f2,cfg_mr2.lfp(iCh)); % Apply filter to data
    %cfg_mr.spkColor = 'w';
end
MultiRaster(cfg_mr2, this_S);
hold on
xline(this_ts)
