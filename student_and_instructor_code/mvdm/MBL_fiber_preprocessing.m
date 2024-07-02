%% Welcome to your first fiber photometry analysis !
% Author: Miriam Janssen
% Date: 07/01/2024
% Citations: Thomas Akam from https://github.com/ThomasAkam/photometry_preprocessing/blob/master/Photometry%20data%20preprocessing.ipynb

% Pre-processing steps: 
% 1. Filtering: to reduce noise and electrical artifacts 
% 2. Detrending: to correct for photobleaching 
% 3. (Omitted) Movement correction: to remove movement artifacts
% 4. Normalization: conversion to dF/F or Z-scoring or both

% -----------------------------------------------------------------------
% manually change file name. your processed data will be saved based on
% this name. 
file_name = 'M514_2024_07-01'; 

% -----------------------------------------------------------------------
% add the vandermeer codebase to your path.
cd('C:\Users\mimia\Documents\Toolboxes\vandermeerlab-replay-da\code-matlab\shared'); % or, wherever your code is located -- NOTE \shared subfolder!
p = genpath(pwd); % create list of all folders from here
addpath(p);

%% load CSV
% folder where data is
cd 'C:\Data\MBL\M514'
data_table = readtable('M514_2024_06_29_Lhemi_pedestal.csv');

%% Save in structure that works for our code base 
% sampling rate is x = samples / second 
sampling_rate = 1/(table2array(data_table(2,1) - data_table(1,1))); 
% 100 Hz or 100 samples/s

% rename variables and store in a time series (ts) struct for manipulation without
% manipulated raw files
FP = ts;
FP.data = table2array(data_table(:,2));
FP.tvec = table2array(data_table(:,1)); 
FP.FS = sampling_rate;

%% Start and End Buffer
% consider adding a buffer
% Photobleaching is exponential and often greatest in the first few
% seconds of recording. This paper recommends removing 2-5 seconds from the
% beginning and end of the recording file 
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7853640/

% Start buffer
% remove first 2 seconds; number of samples to remove = 2/0.0002;
%FP(1:2/0.0002) = []; 
%time(1:2/0.0002) = [];

%% Raw Signal
% This plots the raw data with different three different x-axis timescales.
sessionTitle = 'CW_';
last_time = length(FP.tvec); %the value of the last time point is how many seconds the recording was
timerange1 = 10*sampling_rate; %datapoint range for 10 s
timerange2= 100*sampling_rate; %datapoint range for 100 s 
timerange3= length(FP.tvec); % datapoint range for all data
time_ranges = [timerange1, timerange2, timerange3]; %in seconds 

figure(1)
for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(FP.tvec(t_range), FP.data(t_range), 'Color', [0 0.5 0])
    %title([num2str(time_ranges(t_i)),' samples'], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end

%consider downsampling....

%% Denoised
% Why: to filter out electrical noise greater than 10 Hz 
%       "recording has large electrical noise artifacts, likely due to the
%       high gain amplifiers in the photodetectors picking up signals from
%       nearby mobile phone. The artifacts are very short pulses and can be
%       greatly reduced by running a median filter before the standard low
%       pass filter. 
% Method: Median filter & Low pass filter
% Note: Temporal dynamics of the biosensor are on the level of subseconds
% (X), so filtering out > 10 Hz signals should be ok. 

% Median filter: remove electrical artifacts 
FP_denoised = medfilt1(FP.data);

% check once
figure(2)
plot(FP.tvec,FP.data);
title('Median Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')

% Butterworth Low pass filter 
% Note: Trisch lab used 20 Hz
% fc = 3; % frequency
% [b,a] = butter(2,[1 20]/(FS/2));
% %[b,a] = butter(2,[1 20]/(FS/2)); % 2nd order
% %freqz(b,a,[],FP_data.acq.Fs)
% FP_denoised= filter(b,a,FP_denoised);
% xlim([5 inf])

% Butterworth Low pass filter 
% Note: Trisch lab used 20 Hz
fc = 20; % frequency
[b,a] = butter(2,fc/(FP.FS/2)); % 2nd order
%freqz(b,a,[],FP_data.acq.Fs)
FP_denoised= filter(b,a,FP_denoised);
xlim([1 inf])

figure(3)
plot(FP.tvec,FP.data);
hold on
plot(FP.tvec,FP_denoised)
hold off
title('20 Hz Butterworth Filtered FP Signal over Median Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')
legend('median','butterworth','Location','northeast')
xlim([5 inf])

%%
t = FP.tvec; 
%% Windowed Detrend using Locdetrend
addpath(genpath('C:\Users\mimia\Documents\Toolboxes\vandermeerlab-replay-da'));

FP_detrend_60s = locdetrend(FP_denoised,FP.FS,[60 0.5]);

figure(6)
plot(t,FP_detrend_60s)
title('Detrended and Filtered FP Signal (10s)')
ylabel('Signal (V)')
xlabel('Time (s)')

%% Normalization for Windowed Detrend (locdetrend)
zF_win_60s = (FP_detrend_60s - mean(FP_detrend_60s))./std(FP_detrend_60s);

% filtered, detrended, and normalized signal (using 60s window for detrend)
figure(7)
plot(t,zF_win_60s)
title('Filtered, Detrended, and Normalized FP Signal (60s)')
ylabel('Signal z-scored (V)')
xlabel('Time (s)')

%% Save Variables 
filename = append(file_name, "processed.mat");
data.t = t;
data.z = F_zscored; % detrended and z-scored | previously FP_z ... changed to z.
data.dF = dF; % detrended dF (100.*FP_detrended./F_expfit;) 
data.zdF = zdF; % dF , z-scored

data.detrend_60s = FP_detrend_60s; % locdetrend detrended F using 5s windows
data.z_60s = zF_win_60s; % detrended and z-scored F 

save(filename, '-struct','data')
