% Figure_6_1.m
% Figure_6_1.m produces a series of different input-output
% functions that can be used for firing-rate model neurons. While the input
% is typically thought of as "current" it is more accurate to use input
% "conductance" since presynaptic activity can control the fraction of open
% channels (and hence the conductance) of a postsynaptic cell, but the
% current produced also depends on the postsynaptic membrane potential and
% hence the activity of other presynaptic cells.
%
% The four different model f-I curves are
% 1) power-law
% 2) sigmoid
% 3) threshold linear with saturation
% 4) binary (step function).
%
% This code was used to produce Figure 6.1 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = 0:0.1:20;    % vector of input values

%% First method, is prefix*(inputs)^power.
power = 1.5;
prefix = 1;
% Produce output function with method 1.
outputs1 = prefix*inputs.^power;

%% Second method is sigmoid with a maximum rmax, an input to generate 
%  half-maximum, input_mid, and a term "sigma" that indicates the range
%  of inputs over which a substantial change in output is generated. A
%  small sigma means a steeper maximal slope.
rmax2 = 100;        % maximum rate for output
input_mid2 = 8;      % at this level of input, the output is rmax/2
sigma = 2;          % inverse of steepness.
% Now produce output function with method 2.
outputs2 = rmax2./(1+exp(-(inputs-input_mid2)/sigma));

%% Third method is threshold linear with maximum
rmax3 = 100;            % maximum rate for output
input_mid3 = 8;         % input for half-maximum rate
range = 6;              % range of inputs from 0 to max rate
gradient = rmax3/range; % gradient of f-I curve
y_intercept = rmax3/2 - gradient*input_mid3;
outputs3 = y_intercept + gradient*inputs;
outputs3 = max(outputs3,0);
outputs3 = min(outputs3,rmax3);

%% Fourth method is binary, step function, equivalent to either of methods 
%  2 or 3, in the limit sigma = 0 or range = 0, respectively.
%  The midpoint which is the threshold, is needed, as is the maximum rate
rmax4 = 100;                            % maximum rate
input_mid4 = 8;                         % value for step-change in rate
outputs4 = rmax4*(inputs>input_mid4);  % returns rmax4 or 0 

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

% Plot the 4 panels aligned horizontally
figure(1)
clf
subplot('Position',[0.065 0.22 0.17 0.67])
plot(inputs,outputs1,'k','LineWidth',3)
xlabel('Input')
ylabel('Output')

subplot('Position',[0.315 0.22 0.17 0.67])
plot(inputs,outputs2,'k','LineWidth',3)
xlabel('Input')

subplot('Position',[0.565 0.22 0.17 0.67])
plot(inputs,outputs3,'k','LineWidth',3)
xlabel('Input')

subplot('Position',[0.815 0.22 0.17 0.67])
plot(inputs,outputs4,'k','LineWidth',3)
xlabel('Input')

% Label the panels A-D
annotation('textbox',[0.00 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.25 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.50 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.75 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','D')
