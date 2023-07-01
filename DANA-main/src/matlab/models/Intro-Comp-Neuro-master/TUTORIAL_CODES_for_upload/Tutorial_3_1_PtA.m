% Tutorial_3_1_PtA.m
% Code produces a spike-triggered average for a simple voltage-based model
% neuron receiving white-noise current input.
%
% The functions expandbin and STA are required by this code and available
% at github.com/Primon23
% 
% The code uses the Adaptive Exponential Leaky Integrate and Fire Model,
% a 2-variable model than can reproduce many types of neural behavior
% see Naud, Marcille, Clopath, Gerstner, Biol. Cybern. 2006
%
%   This code is a solution for Tutorial 3.1, Part A of Chapter 3 in the
%   textbook
%   An Introductory Course in Computational Neuroscience
%   by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% List of Cell Parameters

G_L = 8e-9;                % Leak conductance (S)
C = 100e-12;                % Capacitance (F) 
E_L = -60e-3;               % Leak potential (V)
V_Thresh = -50e-3;         % Threshold potential (V)
V_Reset = -80e-3;          % Reset potential (V)
deltaT = 2e-3;             % Threshold shift factor (V)
tauw = 50e-3;             % Adaptation time constant (s)
a = 10e-9;                  % adaptation recovery (S)
b = 0.5e-9;                % adaptation strength (A)

Vmax = 50e-3;              % level of voltage to detect a spike

%% Simulation set-up

dt = 2e-5;               % dt in ms
I0=-0.5e-9;        % Minimum applied current
Irange= 1e-9;     % Range of applied currents  (max - min)
Nsteps = 40000;     % Number of distinct applied currents
Ival = I0 + Irange*rand(1,Nsteps);  % Randomized set of values
steplength = 0.005;         % Durations of constant current
tmax = Nsteps*steplength;   % Total time for simulation
nsteplength = round(steplength/dt); % Number of simulation time-steps when current is held constant
tvector = 0:dt:tmax;                % Vector of time-points
Iapp = zeros(size(tvector));        % Vector to hold the applied current

% The following for loop fills the applied current vector, Iapp, with the
% random set of current values, with the value held constant for a number
% of time bins given by nsteplength
for step = 1:Nsteps;                    % Loop through each value of current
    istart = (step-1)*nsteplength+1;    % Time-step to start that value
    istop = step*nsteplength;           % Time-step to end that value
    Iapp(istart:istop) = Ival(step);    % Constant value for the range of time-points
end

v = zeros(size(tvector));       % initialize voltage
v(1) = E_L;
w = zeros(size(tvector));       % initialize adaptation variable
spikes = zeros(size(tvector));

for j = 1:length(tvector)-1     % simulation of tmax

    if ( v(j) > Vmax )          % if there is a spike
        v(j) = V_Reset;         % reset the voltage
        w(j) = w(j) + b;        % increase the adaptation variable by b
        spikes(j) = 1;          % Record the spike as a "1" in the time-bin
    end
    
    % next line integrates the voltage over time, first part is like LIF
    % second part is an exponential spiking term
    % third part includes adaptation
    v(j+1) = v(j) + dt*( G_L*(E_L-v(j) + deltaT*exp((v(j)-V_Thresh)/deltaT) ) ...
       - w(j) + Iapp(j))/C;

   % next line decys the adaptation toward a steady state in between spikes
    w(j+1) = w(j) + dt*( a*(v(j)-E_L) - w(j) )/tauw;
    
end
            
% Now downsample the stimulus and response to 1ms bins using the online
% function expandbin
newdt = 0.001;                      % new time-bin of 1ms
spikes = expandbin(spikes,dt,newdt);
spikes(find(spikes)) = 1;   % replace fractions with a "1" for a spike
Iapp = expandbin(Iapp,dt,newdt);  % The function expandbin does the downsampling

% Next line calculates the spike-triggered average using the online
% function STA
[sta, tcorr] = STA(Iapp, spikes, newdt);
 
%% Plot the figures
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure()
plot(tcorr,sta)
xlabel('Spike lag (ms)')
ylabel('Stimulus strength')
set(gca,'YTick',[ ])
