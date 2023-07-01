% bistable_percept.m
%
% Firing rate model for two populations with recurrent excitation,
% coupled by cross-inhibition. Units have sigmoidal firing-rate curves.
%
% The simulation is long as multiple transitions between different states
% of activity are simulated and recorded so that the distribution of
% dwelling times in a state can be produced.
%
% This code was used to generate Figure 7.8 of
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

figure(1)
clf

%% Set up the simulation parameters
dt = 0.0005;                % Time step to use
tmax = 10000;                   % Maximum time for simulation
tvec = dt:dt:tmax;          % Vector of time points
Nt = length(tvec);          % Number of time points
Ntrials = 1;              % Number of different simulation trials

Nmax = 10000;

transition_threshold = 40;    % Threshold indicates decision is made

%% Parameters that vary across conditions are given in vectors, with each
%  element of the vector to be used on one corresponding trial.
%  First trial is depression only.
%  Second trial is facilitation only.
%  Third trial has both facilitation and depression.


%% Synaptic response parameters

Wee = 5;          % Strength of synaptic feedback.
Wei = -3;        % Strength of synaptic feedback.
W = [Wee Wei; Wei Wee];

% Time constant for synaptic depression, zero means no depression
tau_adapt = 0.25;

% Single unit parameters for firing-rate model
tau_r = 0.010;          % time constant for changes in rate
rmax = 50;              % maximum rate
Ithresh = 15;           % input needed for half of maximum rate
Isigma = 3;             % range of inputs to significantly change rate

%% Stimulus properties that can vary per trial
stim_amp = 8;                   % Mean stimulus amplitude
stim_diff = 0;                  % Amplitude difference between 2 stimuli
t_on = 1;                       % Time of stimulus onset
stim_length = tmax-t_on;        % Duration of stimulus
t_off = t_on + stim_length;     % Time for stimulus offset

% Because the stimulus noise scales to large values with small dt, it is
% computationally safer to keep the noise timescale fixed at a larger value
% than dt.
noise_time_scale = 0.005;           % 5ms is smaller than most cellular timescales
% noise_bin is the number of time steps in the noise timescale, for which
% the noise input will be identical
noise_bin = round(noise_time_scale/dt);

% N_noise_bins is number of distinct noise values
N_noise_bins = ceil(stim_length/noise_time_scale);

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
% LineStyle_vec is to allow different linestyles within a loop
% Note { } is used for a vector of character strings
LineStyle_vec = {'-', '--', ':'};

% Set the stimuli such that stim_amp is the mean and the difference
% (stim_amp1 - stim_amp2) is obtained from stim_diff_vec
% stim_amp1 is amplitude of stimulus to Unit 1
stim_amp1 = stim_amp + 0.5*stim_diff;
% stim_amp2 is amplitude of stimulus to Unit 2
stim_amp2 = stim_amp - 0.5*stim_diff;

mean_r = zeros(2,Nt);   % To accumulate trial-averaged rate
% To accumulate times to threshold, trial-averaged rates will be
% stored in reverse time order.
mean_reverse_r = zeros(2,Nt);
% To count trials that reach each time-point for later
% normalization to find the mean rates since simulations will stop
% when a decision is made and the total number of trials is no
% longer the correct normalization.
trial_counts = zeros(1,Nt);
state_time = zeros(1,Nmax);   % To store times to threshold.

% The variable "counts" is to record number of times each Unit is the
% "winner".
counts = zeros(1,2);

% Loop through trials of identical stimuli to accumulate statistics
stim = zeros(2, Nt);        % Initialize vector for the stimulus
% then add stimulus to appropriate time points.
stim(1,round(t_on/dt):round(t_off/dt)) = stim_amp1;
stim(2,round(t_on/dt):round(t_off/dt)) = stim_amp2;
stim(1,round(t_on/dt):round((t_on+0.25)/dt)) = 0;

% part 1 is determnistic with adaptation
% part 2 is purely noise-driven transitions without adaptation
% part 3 has a combination of noise and weak adaptation
for part = 1:3
    switch part
        case 1
            I_Noise = 0.0;      % amplitude of noise in input current
            adapt_step = 0.01;  % Strength of adaptation
            Wee = 5;         	% Strength of synaptic self-excitation.
            Wei = -3;        	% Strength of synaptic cross-inhibition.
            W = [Wee Wei; Wei Wee];
        case 2
            I_Noise = 11;       % amplitude of noise in input current
            adapt_step = 0;     % Strength of adaptation
            Wee = 1.5;          % Strength of synaptic self-excitation.
            Wei = -1.5;         % Strength of synaptic cross-inhibition.
            W = [Wee Wei; Wei Wee];
        case 3
            I_Noise = 8;        % amplitude of noise in input current
            adapt_step = 0.007; % Strength of adaptation
            Wee = 5;            % Strength of synaptic self-excitation.
            Wei = -3;           % Strength of synaptic cross-inhibition.
    end
    W = [Wee Wei; Wei Wee];     % set up connectivity matrix
    
    Itot = zeros(2,Nt);
    I_adapt = zeros(2,Nt);
    
    % Below is the formula for the firing rate response function that will
    % be used in both the steady state curves and the simulation.
    % This is a function of "x" where "x" is the value of the
    % variable (or vector of values) sent to the function.
    r_of_I = @(x) rmax./(1+exp(-(x-Ithresh)/Isigma));
    
    r = zeros(2,Nt);        % Initialize vector for firing rates
    
    % Now produce a vector of stimulus-independent noise for each cell
    I_Noise_vec = I_Noise*randn(2,Nt)*sqrt(tau_r/dt);
    
    %% Now simulate the model as a function of time
    % In models with multiple variables, each variable is updated using
    % the values of the set of variables on the previous time step
    
    i = 1;              % i denotes the time-step
    % A while loop is used to exit the trial once a sufficient number of
    % transitions is produced
    transition_number = 0;      % counts number of state transitions
    state = 0;                  % state will switch back and forth between 0 and 1
    while  ( ( i < Nt) && ( transition_number < Nmax) ) % Loop through time points
        i = i+1;            % increment index of time-point
        
        % Total current is:
        % feedback current + stimulus current - adapation current
        Itot(:,i) = W'*r(:,i-1) + stim(:,i) - I_adapt(:,i-1);
        
        % Now update the firing rate using the equation
        % dr/dt = -(r + r(Vss) )/tau_r
        r(:,i) = r(:,i-1)*(1-dt/tau_r) + ...
            dt/tau_r*r_of_I(Itot(:,i)+I_Noise_vec(:,i));
        
        % Update the adaptation current using the equation
        % dI/dt = -I/tau_adapt + r*adapt_step;
        I_adapt(:,i) = I_adapt(:,i-1)*(1-dt/tau_adapt) + ...
            r(:,i-1)*adapt_step;
        
        %% The following section detects and records times of state transitions
        % First look for transitions into state 1
        if ( state ~= 1 )
            if ( r(1,i) > transition_threshold )    % state 1 has high r1
                state = 1;
                if ( transition_number > 0 )        % find duration of last state
                    state_time(transition_number) = ...
                        dt*i - last_transition_time;
                end
                
                % just maintain this transition time for next state
                % duration calculation
                last_transition_time = dt*i;
                transition_number = transition_number + 1;
            end
        end
        % Then detect transitions into state 2
        if ( state ~= 2 )
            if ( r(2,i) > transition_threshold )    % state 2 has high r2
                state = 2;
                if ( transition_number > 0 )        % find duration of last state
                    state_time(transition_number) = ...
                        dt*i - last_transition_time;
                end
                
                % just maintain this transition time for next state
                % duration calculation
                last_transition_time = dt*i
                transition_number = transition_number + 1
            end
        end
        
    end;        % go to next time point
    
    
    %% Now plot the individual trial results
    
    % Given the length of simulations, the time-points are sparsely sampled
    figure(1)
    subplot(2,3,part)
    plot(tvec(25:25:end),r(1,25:25:end),'k');     % plot firing rate as a function of time
    hold on
    plot(tvec(25:25:end),r(2,25:25:end),'Color',[0.5 0.5 0.5]);    % plot firing rate as a function of time
    xlabel('Time (s)')            % Label x-axis on last trial
    ylabel('Rate (Hz)')
    axis([0 10 0 60])         % Set the range of each axis
    set(gca,'Layer','top')          % Ensures axes lines are visible on top of the plot
    
    % Then plot the histogram of state durations
    subplot(2,3,part+3)
    [y x ] = hist(state_time(1:transition_number),31)
    plot(x,y,'k')
    xlabel('Duration (s)')
    ylabel('No. of occurrences')
    drawnow
end

% Finally label the panels on the figure
annotation('textbox',[0 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.34 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.62 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0 0.51 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','D')
annotation('textbox',[0.34 0.51 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','E')
annotation('textbox',[0.62 0.51 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','F')
