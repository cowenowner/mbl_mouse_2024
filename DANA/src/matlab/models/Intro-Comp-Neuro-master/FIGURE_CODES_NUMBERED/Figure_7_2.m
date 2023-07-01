% Figure_7_2.m
%
% Plots dr/dt as a function of r for different values of recurrent
% feedback, W.
%
% This code is used to produce Figure 7.2 in the book:
% An Introductory Course in Computational Neuroscience,
% by Paul Miller.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_in = 0:0.1:100;           % Set of rates for the x-axis

rmax = 100;                 % Maximum rate in f-I curve
I_sigma = 20;               % Range of inputs to change rate (inverse steepness)
Ith = 50;                   % Input for half-maximum rate
% Define the firing rate curve as a sigmoid function
f_of_r = @(x) rmax./(1+(exp(-(x-Ith)/I_sigma)));

Wvals = [1 1.06 1.12 1.18]; % Set of values of feedback, W, to use one at a time
Nplots = length(Wvals);     % Number of subplots, one per value of W

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(2)
clf
figure(1)
clf

%% Loop through different values of W, calculated dr/dt as a function of r
%   then plot the resulting curve for each value of W
for i = 1:Nplots
    Wrec = Wvals(i);            % Value of W for this plot
    
    input = Wrec*r_in;          % Input is recurrent feedback only
    r_out(i,:) = f_of_r(input); % Output rate as a function of input
    
    % In figure(1) plot r_out as a function of r_in and compare to the
    % straight line where r_out = r_in
    figure(1)
    subplot('Position',[0.12+0.3*(i-1) 0.2 0.25 0.7])    
    plot(r_in,r_out(i,:));
    hold on
    plot(r_in,r_in)
    
    % In figure(2) calculate dr/dr from the dynamical equations as:
    % dr/dt = (-r + f(I))/tau
    figure(2)    
    drdt(i,:) = (-r_in + r_out(i,:))/tau;
    
    % Now plot dr/dt versus r
    subplot('Position',[0.1+0.95/Nplots*(i-1) 0.2 2/(3*Nplots) 0.7])
    plot(r_in,drdt(i,:),'k','LineWidth',3)
    hold on
    plot(r_in,zeros(size(r_in)),'k','LineWidth',1)
    xlabel('r (Hz)')
    if ( i == 1 )
        ylabel('dr/dt (Hz/sec)')
    end
    axis([0 rmax -1500 1500])
    title(strcat('W = ',num2str(Wrec)))
    
end
% Label the subplots A-D.
annotation('textbox',[0.01 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.28 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.515 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.75 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','D')

