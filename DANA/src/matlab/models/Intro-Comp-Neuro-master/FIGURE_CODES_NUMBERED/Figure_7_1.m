% Figure_7_1.m
%
% A plot of dr/dt as a function of r for a bistable system produced by a
% firing rate model unit with a sigmoidal firing rate curve and excitatory
% recurrent feedback.
%
% Trajectories of r vs t from different initial conditions are also
% produced.
%
% This code was used to produce Figure 7.1 in the textbook 
%
% An Introductory Course in Computational Neuroscience,
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;                  % Clear all prior variables and parameters from memory

r_in = 0:0.1:100;       % vector of rates to use on x-axis

%% Details of unit's f-I curve (firing rate response to input current)
rmax = 100;             % Maximum rate     
Ith = 50;               % Input for half-maximum rate
I_sigma = 20;           % Range of input to produce change in rate
tau = 0.010;            % Time constant for rate changes

% The f-I curve is produced as an inline function of "x" which would be
% passed to it as the input current.
f_of_r = @(x) rmax./(1+(exp(-(x-Ith)/I_sigma)));

Wrec = 1;               % Recurrent feedback strength

x = Wrec*r_in;          % Input current is linear in firing rate
r_out = f_of_r(x);      % Steady state rate depends on input current

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

% figure(1) can be used to see where the fixed points are, such that
% r_in = r_out.
figure(1)
clf
plot(r_in,r_out);
hold on
plot(r_in,r_in)

% figure(2) becomes Figure 7.1 in the textbook, a plot of dr/dt, which is
% equal to zero when r_in = r_out at the fixed points.
figure(2)
clf
drdt = (-r_in + r_out)/tau;
subplot('Position',[0.12 0.2 0.35 0.7])
plot(r_in,drdt,'k','LineWidth',3)
hold on
plot(r_in,zeros(size(r_in)),'k','LineWidth',1)
xlabel('r (Hz)')
ylabel('F(r) = dr/dt (Hz/sec)')

%% Now simulate from many different initial conditions
tmax = 0.2;             % Maximum time to simulate
dt = 0.001;             % Time-step 
t = 0:dt:tmax;          % Vector of time-points
Nt = length(t);         % Number of time-points

% Set up the plot before looping through initial conditions
subplot('Position',[0.6 0.2 0.35 0.7])
for r0 = 0:4:100;       % Set of initial conditions
    r = zeros(1,Nt);    % Set up firing rate vector
    r(1) = r0;          % Initialize firing rate
    for i = 2:Nt;       % Now loop through time
        rinf = f_of_r(r(i-1)*Wrec); % Steady state firing rate
        % Update rate with the exponential Euler method
        r(i) = rinf + (r(i-1)-rinf)*exp(-dt/tau);
    end
    plot(t,r,'k');      % Plot firing rate versus time (Fig. 7.1B)
    hold on
end
ylabel('r (Hz)')
xlabel('t (sec)')

% Label panels A and B in Fig. 7.1
annotation('textbox',[0.00 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.50 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')

