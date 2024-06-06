%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulate the vHC and evoke responses in the mPFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2022
% Edited by SS and AV 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%%%%%%%%%%%%%%% SINGLE PULSE STIM %%%%%%%%%%%%%%%%%%
%Save Stim file path and name input var
Stim_path='C:\SGL_DATA\vHPC_stim_mPFC_excitability_project\10945';
Stim_fileName='params_from_stim_exp_12_21_23_Rat10947.mat';

%Load Previous Stim File-no need to run the following parameters then
% load('C:\SGL_DATA\vHPC_stim_mPFC_excitability_project\params_from_stim_exp_05_12_23.mat');

%Stimulation Parameters
PRM.biphasic = true; % set to false if using old WPI
PRM.biphasic_both_positive = false; % Both postive ONLY if using old WPI.
PRM.time_between_phases_ms = 0; % need some time - say .05ms if using the old orange WPI. 

PRM.pulse_duration_ms = 0.5;
PRM.n_pulses_per_burst = 1;
PRM.time_between_bursts_sec = 30; %s

PRM.PLOT_IT=false;

%Amplitude Range
PRM.output_current_min=1.23; % This is equal to 100uA based on the translation
PRM.output_current_max=6.83; %This is equal to 600uA 
PRM.output_current_step=0.228; %20uA

% Single Stim Arrays...
PRM.output_current_array_in_order=PRM.output_current_min:PRM.output_current_step:PRM.output_current_max; %Diff amplitudes
PRM.output_current_array=PRM.output_current_array_in_order(randperm(length(PRM.output_current_array_in_order)));%Randomize the stim amplitudes
PRM.n_trials = length(PRM.output_current_array);
pulse_sequence_sec=0.1; %Pulse second
%SAVE STIM FILE
save(fullfile(Stim_path,Stim_fileName))
%%
%%%%%%%%%%%%%%%%%% ONLY RUN FROM HERE %%%%%%%%%%%%%%%%%%
tic;
for iT = 1:PRM.n_trials
    NIDQ_Stimulate_IO(pulse_sequence_sec, 'output_voltage', PRM.output_current_array(iT), ...
        'individual_pulse_duration_sec', PRM.pulse_duration_ms/1000, ...
        'biphasic_pulses', PRM.biphasic, ...
        'biphasic_interval_sec', PRM.time_between_phases_ms/1000, ...
        'biphasic_both_positive', PRM.biphasic_both_positive,...
        'PLOT_IT',PRM.PLOT_IT);
    fprintf('%d - %d\n ',iT,PRM.output_current_array(iT));
    pause(PRM.time_between_bursts_sec)
end
toc;
fprintf('Run End Time: %s\n', datestr(now,'HH:MM:SS.FFF'));

%%
%%%%%%%%%%%%%%%%%% BURST CODE %%%%%%%%%%%%%%%%%%

BURST_PRM.biphasic = true; % set to false if using old WPI
BURST_PRM.biphasic_both_positive = false; % Both postive ONLY if using old WPI.
BURST_PRM.pulse_duration_ms = 0.5;
BURST_PRM.time_between_phases_ms = 0; % need some time - say .05ms if using the old orange WPI. 
BURST_PRM.n_pulses_per_burst = 5;
BURST_PRM.intra_burst_interval_ms = 2.5; % 400 Hz
BURST_PRM.time_between_bursts_sec = 2; % for experiment should be 2s
BURST_PRM.n_trials = 30;
%REPLACE THE VOLTAGE WITH 70% of MAX VOLTAGE

BURST_PRM.output_voltage = 3.93; % This is equal to 300uA based on the translation - Abhi
%BURST_PRM.output_voltage = 5.76; % This is equal to 500uA based on the translation - Abhi
% SIMPLE STIM EXPERIMNET...
pulse_sequence_sec = 0:BURST_PRM.intra_burst_interval_ms/1000:(1 + BURST_PRM.n_pulses_per_burst*BURST_PRM.intra_burst_interval_ms/1000);
pulse_sequence_sec = pulse_sequence_sec(1:BURST_PRM.n_pulses_per_burst);

%%
%%%%%%%%%%%%%%%%%% ONLY RUN FROM HERE %%%%%%%%%%%%%%%%%%
tic;
for iT = 1:BURST_PRM.n_trials
    NIDQ_Stimulate_IO(pulse_sequence_sec, 'output_voltage', BURST_PRM.output_voltage, ...
        'individual_pulse_duration_sec', BURST_PRM.pulse_duration_ms/1000, ...
        'biphasic_pulses', BURST_PRM.biphasic, ...
        'biphasic_interval_sec', BURST_PRM.time_between_phases_ms/1000, ...
        'biphasic_both_positive', BURST_PRM.biphasic_both_positive);
    fprintf('%d ',iT);
    pause(BURST_PRM.time_between_bursts_sec)
end



pause (330);
for iT = 1:BURST_PRM.n_trials
    NIDQ_Stimulate_IO(pulse_sequence_sec, 'output_voltage', BURST_PRM.output_voltage, ...
        'individual_pulse_duration_sec', BURST_PRM.pulse_duration_ms/1000, ...
        'biphasic_pulses', BURST_PRM.biphasic, ...
        'biphasic_interval_sec', BURST_PRM.time_between_phases_ms/1000, ...
        'biphasic_both_positive', BURST_PRM.biphasic_both_positive);
    fprintf('%d ',iT);
    pause(BURST_PRM.time_between_bursts_sec)
end
toc;