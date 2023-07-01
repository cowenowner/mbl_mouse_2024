function OUT = Looking_at_burst_activity_using_Burst_Detector_fucntion
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DIRS
OUT = [];
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();

% SES = LK_Session_Info();
max_ISI_ms = 10;
min_spikes_per_burst = 3;

[first_in_burst, durations, n_per_burst, burst_id] = Burst_detector_Abhi(TS, max_ISI_ms*1000, min_spikes_per_burst);