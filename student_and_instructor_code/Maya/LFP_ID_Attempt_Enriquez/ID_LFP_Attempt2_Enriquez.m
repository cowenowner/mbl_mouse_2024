%% Identifying LFP Waveforms
% This script is the second attempt to make a filter that will find sharp-wave
% ripples in randomly generated waveform data.

%% Creating random dataset
fake_data = randn(20000,1);
plot (fake_data);

%% Creating a time vector
Fs = 1000 %Defining sampling rate
tvec = 0:1:length(fake_data)-1; %Making an empty time vector same length as data
tvec = tvec *(1/Fs) %Converting sampling rate into seconds (for sanity)
plot(tvec)%To confirm it actually converted to seconds
plot(tvec, fake_data)%To confirm it actually converted to seconds, again

%% Normalizing the data using z-score so we can generalize across trials

z = zscore(fake_data)

plot (tvec, zscore)
%% Creating a filter
% With help from MVDM and Abhi!

%Filter 1, in Radians (an example; pls don't use this)
%fs = 1000;
%fcutlow = 140;
%fcuthigh = 200;
%[b,a] = butter(10,[fcutlow, fcuthigh]/(fs/2));
%Filt_filt_signal = filtfilt(b, a, fake_data);
%hold on
%plot(Filt_filt_signal,'r')


%Filter 2 in Hz (much nicer)
d_ripple = designfilt('bandpassiir','FilterOrder',12, ...
                    'HalfPowerFrequency1',140,'HalfPowerFrequency2',200, ...
                    'SampleRate',1000); %Adjust for actual sampling rate

filt_data = filtfilt(d_ripple,z);
plot(filt_data);

%% Creating an envelope
% With help from MVDM!

env = abs(hilbert(filt_data));
plot(env,'g')

%% Finding peaks of the envelope (potential sharp wave ripples)
[pks,locs, w, p] = findpeaks(env);% locs is locations, pks is local peak

keep_idx_w = w <15 %Make an index where width must be greater than 15 
w(1:10) %Check first 10 scores of w
keep_idx_w(1:10) %Check the logical array of w; should be 1 for those above 15, 0 for those below

pks = pks(keep_idx); locs = locs(keep_idx) %Synchronize data with peaks of width greater than 15 


