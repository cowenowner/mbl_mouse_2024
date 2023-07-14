%% Odor encoding
% Question: Is odor identity encoded in dorsal CA1 neuron activity

%% set the path
% Command summary goes here
restoredefaultpath;
addpath(genpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\vandermeerlab\code-matlab\shared'));
addpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\mvdm');
addpath('C:\Users\Administrator\Documents\GitHub\mbl_mouse_2023\student_and_instructor_code\Maya\UtilityFunctions');


%% load the data
cd('E:\Data\HC-recordings\SecondRecordingWithAutomatedOdors\preprocessed');
load AllSpikes.mat
load odor_events.mat

% select only those neurons that are actually in dorsal CA1
S = ConvertSPtoS(SP);

%% optioanl viz
cfg_mr = [];
cfg_mr.evt = evt;
MultiRaster(cfg_mr, S);

%% Timing of stimulus

% set start and end times for each odor -- requires "windowsize" parameter
windowsize = 

% for each odor (1-3)

event_time1 = evt.t{1}; 
event_time2 = evt.t{2};
event_time3 = evt.t{3};

% loop through each trial for that odor (start and end time) and get vector of spike counts

% arrange each vector in its place in matrix of neurons x trials

% plot the matrix (imagesc)

%% 

% for each neuron,

% plot firing rate as a function of time-to-event (peri-event time
% histogram)