% Figure_5_2.m
% Plots the time-course of a synapse with separate rise and decay times.
%
% This code is used to produce Figure 5.2 of the textbook,
% "An Introductory Course in Computational Neuroscience"
% by Paul Miller, Brandeis University, 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.0001;            % Time step
t = 0:dt:0.3;           % Time vector
gsyn = zeros(size(t));  % Vector to store total conductance
tspikes = [0.050 0.100 0.150];  % Times of spikes

tau_rise = 0.005;       % Synaptic rise time constant
tau_decay = 0.050;      % Synaptic decay time constant
delta_g = 1e-9;         % Peak change in conductance
% The parameter K below is a normalization to ensure the peak change from a
% spike is delta_g
K = (tau_decay/(tau_rise+tau_decay))*...
    (tau_rise/(tau_rise+tau_decay))^(tau_rise/tau_decay);

% Loop through the spikes
for n = 1:length(tspikes)
    ispike = round(tspikes(n)/dt);      % Time-point of the spike
    % Now update the conductance for all times after the spike using 
    % Eq. 4.4
    gsyn(ispike+1:end) = gsyn(ispike+1:end) + (delta_g/K)*...
        exp(-(t(ispike+1:end)-t(ispike))/tau_decay).* ...
        (1 - exp(-(t(ispike+1:end)-t(ispike))/tau_rise) );
end
%% Now plot the figure
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf
plot(t*1000,gsyn*1e9,'k')
xlabel('Time (msec)')
ylabel('G_{syn} (nS)')
