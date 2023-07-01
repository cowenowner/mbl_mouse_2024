% spatial_AELIF_randomI.m
%
% This code is used to produce a neuron with a spatio-temporal receptive 
% field, by providing inputs in a single spatial dimension weighted by a
% Gabor filter (a cosine multiplied by a Gaussian) to an 
% Adaptive Exponential Leaky Integrate and Fire Model neuron 
% (a 2-variable model than can reproduce many types of neural behavior
% see Naud, Marcille, Clopath, Gerstner, Biol. Cybern. 2006)
%
% This code is used to produce Figure 3.5 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% List of Cell Parameters

G_L = 8e-9;                 % Leak conductance (S)
C = 100e-12;                % Capacitance (F) 
E_L = -60e-3;               % Leak potential (V)
V_Thresh = -50e-3;          % Threshold potential (V)
V_Reset = -80e-3;           % Reset potential (V)
deltaT = 2e-3;              % Threshold shift factor (V)
tauw = 50e-3;               % Adaptation time constant (s)
a = 40e-9;                  % adaptation recovery (S)
b = 1e-9;                   % adaptation strength (A)

Vmax = 50e-3;               % level of voltage to detect a spike

%% Simulation set-up

dt = 2e-5;                  % dt in ms
I0=-0.5e-9;                 % Minimum applied current
Irange= 1e-9;               % Range of applied currents  (max - min)
Nsteps = 40000;             % Number of time-steps of distinct applied currents
Nspatial = 40;              % Number of spatially distinct values per time-step
Ival = I0 + Irange*rand(Nspatial,Nsteps);  % Randomized set of values
steplength = 0.005;         % Durations of constant current
nsteplength = round(steplength/dt); % No. of time points per current step
tmax = Nsteps*steplength;   % Total time for simulation
tvector = 0:dt:tmax;        % Vector of time points
Iapp = zeros(size(tvector));    % Initialize input current vector

%% Set up spatial input vector
x = -(Nspatial-1)/2:(Nspatial-1)/2;     % Vector of spatial x-values
% The input weights below are a Gabor function with a decaying spatial
% oscillation
input_weights = cos(4*pi*x/Nspatial) ...
    .*exp(-4*x.*x/(0.25*Nspatial*Nspatial));

for step = 1:Nsteps;                    % loop through steps of fixed current
    istart = (step-1)*nsteplength+1;    % begin step of fixed current
    istop = step*nsteplength;           % end step of fixed current
    Iapp(istart:istop) = input_weights*Ival(:,step);  % value of fixed current
end

v = zeros(size(tvector));       % initialize voltage
v(1) = E_L;                     % begin at leak potential
w = zeros(size(tvector));       % initialize adaptation variable as zero
spikes = zeros(size(tvector));  % initialize spike vector as zero (no spikes)

for j = 1:length(tvector)-1     % simulation through all time points

    if ( v(j) > Vmax )          % if there is a spike
        v(j) = V_Reset;         % reset the voltage
        w(j) = w(j) + b;        % increase the adaptation variable by b
        spikes(j) = 1;          % record the spike time
    end
    
    % next line integrates the voltage over time, first part is like LIF
    % second part is an exponential spiking term
    % third part includes adaptation
    v(j+1) = v(j) + dt*( G_L*(E_L-v(j) + deltaT*exp((v(j)-V_Thresh)/deltaT) ) ...
       - w(j) + Iapp(j))/C;

   % next line decys the adaptation toward a steady state in between spikes
    w(j+1) = w(j) + dt*( a*(v(j)-E_L) - w(j) )/tauw;
    
end

%% Now downsample the data

newdt = 0.001;                          % Use 1ms bins
spikes = expandbin(spikes,dt,newdt);    % New spike vector with 1ms bins
newNt = length(spikes);                 % New length of spike vector
newIapp = zeros(Nspatial,newNt);        % Define a new input vector
newnsteplength = round(steplength/newdt); % New number of 1ms bins per step

for step = 1:Nsteps;                    % Loop through steps of constant current
    istart = (step-1)*newnsteplength+1; % first time point with 1ms steps
    istop = step*newnsteplength;        % last time point with 1ms steps
    
    % generate the input vector as before, but with bins of 1ms, not of dt
    newIapp(:,istart:istop) = Ival(:,step)*ones(1,newnsteplength);
end

% calculate the spike-triggered average using the online function
% STA_spatial.
[sta, tcorr] = STA_spatial(newIapp, spikes, newdt); 

%% Now plot the spike-triggered average
% Default plotting parameters 
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure()
imagesc(fliplr(sta));       % reverses time-axis to plot STA
colormap(gray)              % grayscale
set(gca,'XTick',[1, 26, 51, 76, 101])
set(gca,'XTickLabel',{'-25' '0', '25', '50', '75'})
xlabel('Spike lag (ms)')
ylabel('Stimulus coordinate')
set(gca,'YTick',[ ])

   