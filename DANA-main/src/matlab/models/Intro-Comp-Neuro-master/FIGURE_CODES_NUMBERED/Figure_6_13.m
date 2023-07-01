% Figure_6_13.m
%
% Firing rate model for one inhibitory population with constant excitation
% and delayed feedback inhibition to produce oscillations.
%
% This code produces Figure 6.13 in the text book:
% An Introductory Course in Computational Neuroscience
% by Paul Miller.
%
clear;

%% Set up the simulation parameters
dt = 0.0001;                 % Time step to use
tmax = 1.5;                   % Maximum time for simulation
tvec = dt:dt:tmax;           % Vector of time points
Nt = length(tvec);           % Number of time points

%% Parameters for firing-rate response curve which is the model by
%  FS Chance and LF Abbott, Prog Brain Res 149:147-155, 2005
%  The model approximates well the output of a leaky integrate-and-fire
%  (LIF) neuron given its steady state voltage. So LIF parameters are
%  needed.
delta_V = 0.030;            % Threshold voltage minus reset voltage
% The parameters a and sigma_V both set the curvature of firing rate curve.
% sigma_V should be the standard deviation of voltage noise and then gets
% scaled (inversely) by a. Here sigma_V does not actually add noise
% fluctuations unless the flag V_noise_on_flag is set to 1.
a = 1;
sigma_V = 0.001;
Vth = -0.050;               % Threshold voltage (V)
E_L = -0.070;               % Leak voltage (V)
G_L = 0.05;                 % Leak conductance (nS)
E_I = -0.065;               % Inhibitory synapse reversal potential (V)
tau = 0.005;                % Time constant for firing-rate to change

%% Parameters that vary across conditions are given in vectors, with each
%  element of the vector to be used on one corresponding trial.
%  First trial is depression only.
%  Second trial is facilitation only.
%  Third trial has both facilitation and depression.


%% Synaptic response parameters

alpha0 = 1;                         % Baseline postsynaptic amplitude.
Wii = 10;                           % Strength of I-to-I synaptic feedback.

tau_sI = 0.005;                     % Synaptic time constant for I-feedback
% Delay for rate change to affect synaptic input
tau_s_delay = 0.003;                % delay in seconds 
i_delay = round(tau_s_delay/dt);    % delay in time-points

%% Stimulus properties that vary per trial

t_on = dt;                          % Time of stimulus onset
stim_length = tmax-t_on;            % Duration of stimulus

% G_input is amplitude of excitatory input to Unit
G_input = 0.15;

G_stim = zeros(1, Nt);        % Initialize vector for the stimulus
% then add stimulus to appropriate time points.
t_off = t_on + stim_length;         % Time for stimulus offset
G_stim(round(t_on/dt):round(t_off/dt)) = G_input;

% Below is the formula for the firing rate response function that will
% be used in both the steady state curves and the simulation.
% This is a function of "x" where "x" is the value of the
% variable (or vector of values) sent to the function.
r_of_Vss = @(x) (x - Vth)./(tau*delta_V*(1-exp(-a*(x-Vth)/sigma_V)));

% To prevent division by zero when x = Vth, use the limit
% value:
r_of_Vss_limit = 1./(tau*delta_V*a/sigma_V);

r = zeros(1,Nt);        % Initialize vector for firing rates
sI = zeros(1,Nt);       % Initialize vector for synaptic inputs

%% Now carry simulate the same model as a function of time
% In models with multiple variables, each variable is updated using
% the values of the set of variables on the previous time step

% A while loop is used to exit the trial once threshold is
% reached
for i =2:Nt
    
    % Calculate the steady-state membrane potential for a
    % leaky-integrate-and-fire cell based on the excitatory and
    % inhibitory synaptic inputs.
    % Vss is a column vector with one value for each cell.
    % The synaptic inputs are obtained by matrix multiplication
    % using the connectivity matrices Wee_mat and Wei_mat, plus
    % the stimulus input conductance in "stim".
    Vss = (E_L*G_L +E_I*Wii*sI(i-1) )./ ...
        (G_L  + Wii*sI(i-1) +G_stim(i));
    
    % Now update the firing rate using the equation
    % dr/dt = -(r + r(Vss) )/tau
    if ( Vss == Vth )
        r(i) = r(i-1)*(1-dt/tau) + ...
            dt/tau*r_of_Vss_limit;
    else
        r(i) = r(i-1)*(1-dt/tau) + ...
            dt/tau*r_of_Vss(Vss);
    end
    
    % Next update the synaptic response variable,s, which depends on
    % the firing rate a few timepoints earlier, using the equation:
    % ds(t)/dt = -s(t)/tau_s + alpha.r(t').( 1 - s(t) )
    % where t' is the earlier time, t' = t - t_delay
    if ( i > i_delay)   % delay before updating the synaptic variable
        sI(i) = sI(i-1)*(1-dt/tau_sI) + ...
            dt*alpha0.*r(i-i_delay).*(1-sI(i-1));
    end
end;        % go to next time point

%% Set up the plotting parameters and plot the results
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

% Plot the results with a time offset because of the initial transient
figure(1)
tmin = 0.1;                 % Initial offset after the transient has ended
plot(tvec-tmin,r,'k')
hold on
axis([0 0.1 0 20])
drawnow
xlabel('Time (sec)')
ylabel('Firing Rate (Hz)')


