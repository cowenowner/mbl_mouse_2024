%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulate the vHC and evoke responses in the mPFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2022
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRM.biphasic = true; % set to false if using old WPI
PRM.biphasic_both_positive = false; % Both postive ONLY if using old WPI.
PRM.time_between_phases_ms = 0; % need some time - say .05ms if using the old orange WPI. 

PRM.pulse_duration_ms = 0.5;
PRM.n_pulses_per_burst = 1;
PRM.time_between_bursts_sec = 0.5; %30s

PRM.output_current_min=1.18; % This is equal to 100uA based on the translation
PRM.output_current_max=4.61; %This is equal to 400uA 
PRM.output_current_step=0.1;

PRM.PLOT_IT=true;
%%%%%%%% SIMPLE STIM EXPERIMNET...
PRM.output_current_array=PRM.output_current_min:PRM.output_current_step:PRM.output_current_max;
PRM.n_trials = length(PRM.output_current_array);

pulse_sequence_sec=0.1; %Pulse second
tic;

% NOTE: we could just make a single train with bursts and send that out.
% This, however, makes it easier to break the train if things go wrong. It
% won't have the long-pulse interval precision though, but that does not
% matter.
for iT = 1:PRM.n_trials
    NIDQ_Stimulate_IO(pulse_sequence_sec, 'output_voltage', PRM.output_current_array(iT), ...
        'individual_pulse_duration_sec', PRM.pulse_duration_ms/1000, ...
        'biphasic_pulses', PRM.biphasic, ...
        'biphasic_interval_sec', PRM.time_between_phases_ms/1000, ...
        'biphasic_both_positive', PRM.biphasic_both_positive,...
        'PLOT_IT',PRM.PLOT_IT);
    fprintf('%d - %d ',iT,PRM.output_current_array(iT));
    pause(PRM.time_between_bursts_sec)
end


toc;
