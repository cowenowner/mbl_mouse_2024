% Figure_1_6.m
%
% Example of constant acceleration to produce a velocity of a.t and a
% position of a.t^2.
%
% This code first produces the analytic result then solves the differential 
% equations by accumulating as in Table 1.2.
%
% This code was used to produce Figure 1.5 in the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First the analytic results, i.e. use the formulae.
dt = 0.01;          % Time-step for plot (in sec)
tmax = 5;           % Maximum time for plot
t = 0:dt:tmax;      % Vector of time-points
a = -10;            % Constant acceleration (gravity)
v = a*t;           % Linearly varying velocity
y = 0.5*a*t.*t;    % Quadratically varying position

%% Next do the simulation as in Table 1.2 but with small time-steps
v_sim = zeros(size(t));
y_sim = zeros(size(t));
for i = 2:length(t)
    v_sim(i) = v_sim(i-1) + dt*a;
    y_sim(i) = y_sim(i-1) + dt*0.5*(v_sim(i-1)+v_sim(i));
end
%% Set up the plotting parameters and plot the results
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf
subplot(2,1,1)
plot(t,v,'k')
hold on
xlabel('Time (sec)')
ylabel('Velocity (ms^{-1})')
plot([0 1 2 3 4],[0 -10 -20 -30 -40],'xk')
axis([0 5 -50 0])
subplot(2,1,2)
plot(t,y,'k')
hold on
xlabel('Time (sec)')
ylabel('Height (m)')
plot([0 1 2 3 4],[0 -5 -20 -45 -80],'xk')
axis([0 5 -130 0])
annotation('textbox',[0.00 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.00 0.49 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')

figure(2)
subplot(2,1,1)
plot(t,v_sim)
hold on
plot(t,v)

subplot(2,1,2)
plot(t,y_sim)
hold on
plot(t,y)
