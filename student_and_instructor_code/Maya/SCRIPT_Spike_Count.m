%% Spike Count
% This script is used to count the number of spikes in a given time period

%% Loading the code
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\mvdm'));
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\Maya'));
addpath(genpath('C:\data\02_7_7_23\preprocessed'));

%% Load the data
datadir = ('C:\data\02_7_7_23\preprocessed'); %Create data directory

load(fullfile(datadir,'AllSpikes.mat')); % Load in experiment data
load(fullfile(datadir,'odor_events.mat')); %Load in odor presentation times
load(fullfile(datadir,'LFPs.mat')); %Load in LFP data HUUUGE

%% Create event time vectors
%Create a vector containing stimulus presentation times
event_time1 = evt.t{1}; 
event_time2 = evt.t{2};
event_time3 = evt.t{3};

%% Determine total spikes during baseline
baseline1 = countspk(0, 600, SP);
exp = countspk(601, 1853, SP);
baseline2 = countspk(1853,2489,SP)

plot(baseline1)
hold on
plot(baseline2)
plot(exp)