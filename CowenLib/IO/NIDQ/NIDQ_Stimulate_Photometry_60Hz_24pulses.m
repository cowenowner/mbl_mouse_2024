%%
%%%%%%%%%%%%%%%%%% BURST CODE %%%%%%%%%%%%%%%%%%

BURST_PRM.biphasic = true; % set to false if using old WPI
BURST_PRM.biphasic_both_positive = false; % Both postive ONLY if using old WPI.
BURST_PRM.pulse_duration_ms = 2;
BURST_PRM.time_between_phases_ms = 0; % need some time - say .05ms if using the old orange WPI. 
BURST_PRM.n_pulses_per_burst = 60;
BURST_PRM.intra_burst_interval_ms = 16.66; % 60 Hz
BURST_PRM.time_between_bursts_sec = 60; % for experiment should be 2s
BURST_PRM.n_trials = 1;
%REPLACE THE VOLTAGE WITH 70% of MAX VOLTAGE

BURST_PRM.output_voltage = 3.48; % This is equal to 300uA based on the translation - Abhi
%BURST_PRM.output_voltage = 5.76; % This is equal to 500uA based on the translation - Abhi
% SIMPLE STIM EXPERIMNET...
pulse_sequence_sec = 0:BURST_PRM.intra_burst_interval_ms/1000:(1 + BURST_PRM.n_pulses_per_burst*BURST_PRM.intra_burst_interval_ms/1000);
pulse_sequence_sec = pulse_sequence_sec(1:BURST_PRM.n_pulses_per_burst);

%%
%%%%%%%%%%%%%%%%%% ONLY RUN FROM HERE %%%%%%%%%%%%%%%%%%

for iT = 1:BURST_PRM.n_trials
    NIDQ_Stimulate_IO(pulse_sequence_sec, 'output_voltage', BURST_PRM.output_voltage, ...
        'individual_pulse_duration_sec', BURST_PRM.pulse_duration_ms/1000, ...
        'biphasic_pulses', BURST_PRM.biphasic, ...
        'biphasic_interval_sec', BURST_PRM.time_between_phases_ms/1000, ...
        'biphasic_both_positive', BURST_PRM.biphasic_both_positive);
    fprintf('%d ',iT);
    pause(BURST_PRM.time_between_bursts_sec)
end
