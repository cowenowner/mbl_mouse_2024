%% Identifying LFP Waveforms
% This script is the first attempt to make a filter that will find sharp-wave
% ripples in randomly generated waveform data. 

%Learning randn() documentation
doc randn

%Learning filtfilt() documentation
doc filtfilt

%% Generating fake dataset
fake_data = randn(10000,1);
plot (fake_data)

%% Attempt A
%Specifying constraints of the filter based on this link: 
% https://stackoverflow.com/questions/62209954/phase-difference-removal-by-using-filtering-butterworth-filtfilt-command%
fs = 1000;
fcutlow = 140;
fcuthigh = 200;
[b,a] = butter(10,[fcutlow, fcuthigh]/(fs/2));
%Butterworth_bandpass_filter  = filter(b,a,fake_data);
Filt_filt_signal = filtfilt(b, a, fake_data);
hold on
plot(Filt_filt_signal,'r')

env = abs(hilbert(Filt_filt_signal));
plot(env,'g')
%% Attempt B
%Specifying constraints of filter based on example in filtfilt doc, Matlab
%Help from Abhi!


d = designfilt("lowpassfir", ...
    'PassbandFrequency',0.18,'StopbandFrequency',0.2, ...
     'PassbandRipple',1,'StopbandAttenuation',60);


d14 = designfilt('lowpassiir','FilterOrder',14, ...
        'StopbandFrequency',0.2,...
         'StopbandAttenuation',60);

d_ripple = designfilt('bandpassiir','FilterOrder',12, ...
                    'HalfPowerFrequency1',140,'HalfPowerFrequency2',200, ...
                    'SampleRate',1000);

filt_fake = filtfilt(d_ripple,fake_data);
%% Attempt C
%Specifying constraints of filter based on example in bandpass doc, Matlab
bandpass(fake_data,[180 200],2500)
doc filtfilt


