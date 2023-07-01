%exponential.m
%
% A plot of membrane potential (y1 or y2) against time for a neuron with a
% leak potential of -70mV and no input current, commencing from two
% different initial values.
%
% This code was used to produce Figure 2.2 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University, 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

tau = 0.010;        % Membrane time constant
t = 0:0.001:0.05;   % vector of time points

E_L = -0.070;       % Leak potential (steady state or resting potential)

y1 = E_L + (-0.050 - E_L)*exp(-t/tau);  % Decay from -50mV 
y2 = E_L + (-0.080 - E_L)*exp(-t/tau);  % Decay from -80mV

% Now clear the current figure then plot the two curves.
% COWEN - now do this the non-analytical/computational way.
clf
plot(t,y1,'k');
hold on
plot(t,y2,'k');
plot(t,E_L*ones(size(t)),'--k') % Plot leak potential as a dashed line
xlabel('Time (sec)')
ylabel('V_m (V)')

% remember, you need tau - that will have to be calculated or optimized
% it is easy to screw up the sign.
y = zeros(size(t));
y(1) = -0.090;
dt = 0.001;
for ii = 2:length(t)
    y(ii) = y(ii-1) + dt*(E_L - y(ii-1))/tau;
end
figure
plot(t,y)
hold on
plot(t,E_L*ones(size(t)),'--k') % Plot leak potential as a dashed line
