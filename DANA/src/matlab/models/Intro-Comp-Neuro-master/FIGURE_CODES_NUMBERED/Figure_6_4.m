% Figure_6_4.m
% Simple firing-rate model of bistability.
% Uses sigmoidal f-I curve so that firing rate saturates.
% Plots firing rate and feedback curves to show fixed points then simulates
% the model with inputs.
% Three cases with different feedback strengths are depicted, corresponding
% to the three regimes of single stable quiescent state; two stable states;
% and single stable high-activity state.
%
% This code was used to produce Figure 6.4 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%% Parameters for sigmoidal firing-rate response curve
rmax = 100;                 % Maximum rate
x_sigma = 0.1;              % Inverse of steepness of response
x_mid = 0.3;                % Default input for 1/2-max response (can vary)
tau = 0.010;                % Time constant for firing-rate to change

%% Parameters that vary across conditions are given in vectors, with each 
%  element of the vector to be used on one corresponding trial.


%% Synaptic response parameters
Wee_vec = [0.2 0.6 1.2];        % 3 different strengths of synaptic feedback.

tau_s = 0.002;                  % Synaptic time constant for AMPA receptors

%% Now set up vectors to generate the firing rate response to a range of 
%  inputs and the dependence of synaptic feedback on a range of rates.
S_range = [0:0.001:1.5];          % Range of synaptic inputs to use    
r_range = [0:0.01:rmax];        % Range of presynaptic rates to use.


%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf
%% Set up the simulation parameters 
dt = 0.0005;                % Time step to use
tmax = 1;                   % Maximum time for simulation
tvec = dt:dt:tmax;          % Vector of time points
Nt = length(tvec);          % Number of time points
Ntrials = 3;                % Number of different simulation trials

r = zeros(Ntrials,Nt);      % Initialize vector for firing rates
s = zeros(Ntrials,Nt);      % Initialize vector for synaptic inputs

t_on = 0.5;                   % Time for stimulus onset
t_off = 0.55;                % Time for stimulus offset
stim_amp = 0.2;             % Amplitude of stimulus
stim = zeros(1, Nt);        % Initialize vector for the stimulus
% then add stimulus to appropriate time points.
stim(round(t_on/dt):round(t_off/dt)) = stim_amp;    

% Loop through trials with different parameters on each trial
for trial = 1:Ntrials;         
    Wee = Wee_vec(trial);       % Value of total synaptic strength
    
    % Below is the formula for the firing rate response function that will 
    % be used in both the steady state curves and the simulation.
    % This is the sigmoid function of "x" where "x" is the value of the
    % variable (or vector of values) sent to the function.
    % This is written as an inline function as it will be used many times.
    r_of_s = @(x) rmax./(1 + exp(-(x-x_mid)/x_sigma));

    % Generate firing rate as a function of total synaptic input
    r_out = r_of_s(S_range);
    
    % Steady state formula for the synaptic response fraction depends 
    % linearly on r in this approximation.
    s_of_r = r_range/rmax;
    
    % Total synaptic input is the synaptic response fraction scaled by
    % total synaptic strength, Wee via S(r) = Wee.s(r).
    S_in = Wee*s_of_r;
        
    % Set up the figure to plot S(r) and r(S).
    figure(1)    
    subplot('Position',[0.12 1.04-0.31*trial 0.33 0.22]) % Position of panel
    plot(S_range,r_out,'k');        % Plot r(S) (S varied r from formula)
    hold on;                        % Do not overwrite
    plot(S_in,r_range,'k--');       % Plot S(r) (r varied S from formula)
    axis([0 1.5 0 rmax])               % Set the range of each axis

    if ( trial == Ntrials )          % If it is not the last trial
        xlabel('Synaptic Input')    % Otherwise label the x-axis
    end
    ylabel('Rate (Hz)')             % Label the y-axis
    if ( trial == 1 )
        legend('rate', 'feedback')
    end

    %% Now carry simulate the same model as a function of time
    % In models with multiple variables, each variable is updated using
    % the values of the set of variables on the previous time step
    for i = 2:Nt;                               % Loop through time points

        % First update the firing rate using the equation
        % dr/dt = -(r + r(S + stim) )/tau 
        % where S = Wee.s and the value of s on the prior time-step is used
        r(trial,i) = r(trial,i-1)*(1-dt/tau) + ...
            dt/tau*r_of_s(Wee*s(trial,i-1)+stim(i));
        
        
        % Next update the synaptic response variable,s, which depends on
        % rate, using the equation
        % ds/dt = (-s + r/rmax)/tau_s
        s(trial,i) = s(trial,i-1)*(1-dt/tau_s) + ...
            dt/tau_s*r(trial,i-1)/rmax;
        
    end;        % go to next time point
    
    %% Now plot the results
    
    % Set up the position of the panel for plotting simulation results
    subplot('Position',[0.62 1.04-0.31*trial 0.36 0.22])
    
    % produce a solid gray rectangle across the region of the stimulus
    fill([t_on t_on t_off t_off],[0 rmax rmax 0],[0.75, 0.75, 0.75], ...
        'EdgeColor','none')
    hold on;        % do not overwrite the gray rectangle
    plot(tvec,r(trial,:),'k');  % plot firing rate as a function of time
    if ( trial == Ntrials )          % if not the last trial
        xlabel('Time (sec)')        % Label x-axis on last trial
    end
    ylabel('Rate (Hz)')
    axis([0 tmax 0 rmax])            % Set the range of each axis
    
end;            % Go to next parameter set and produce a new trial


% Finally label the subfigures
annotation('textbox',[0.00 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A1')
annotation('textbox',[0.5 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A2')
annotation('textbox',[0.00 0.68 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B1')
annotation('textbox',[0.5 0.68 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B2')
annotation('textbox',[0.00 0.38 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C1')
annotation('textbox',[0.5 0.38 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C2')
