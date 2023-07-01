% integrator_cartoon.m
%
% Code demonstrating how a perfect integrator produces parameteric working
% memory.
%
% No simulations carried out -- "cumsum" is used to produce integration
%
% This code is used to produce Figure 6.7 in the text book
% An Introductory Course in Computational Neuroscience 
% by Paul Miller.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.001;             % time-step
tmax = 3;               % maximum time to plot
tvec = 0:dt:tmax;       % vector of time points
    
stimsvec = [10 20 30];          % set of differnet stimulus strengths
Nstims = length(stimsvec);      % no. of different stimuli 
stim_on = 1;                    % time to commence stimulus
stim_off = 1.5;                 % time to end stimulus
n_on = 1+round(stim_on/dt);     % time-point number to start stimulus
n_off = 1+round(stim_off/dt);   % time-point number to end stimulus

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
% LineStyle_vec is to allow different linestyles within a loop
% Note { } is used for a vector of character strings
LineStyle_vec = {'-', '--', ':'};

figure(1)
clf

%% Now loop through with each different stimulus strength
for i = 1:Nstims
    % first produce the stimulus vector and plot it
    stim = zeros(size(tvec));           
    stim(n_on:n_off) = stimsvec(i);
    subplot('Position',[0.16 0.59 0.8 0.35])
    plot(tvec,stim,'k','linestyle',LineStyle_vec{i})
    hold on
    ylabel('Input Strength')
    axis([0 tmax -5 35])
    
    % next use "cumsum" to calculate the cumulative sum of the stimulus,
    % which given the small dt is equivalent to its integral so a plot of
    % cumsum(stim)*dt is a plot of the integrated stimulus
    subplot('Position',[0.16 0.11 0.8 0.35])
    plot(tvec,dt*cumsum(stim),'k','linestyle',LineStyle_vec{i})
    hold on
    ylabel('Integrated Output')
    xlabel('Time (sec)')
    axis([0 tmax -2.5 17.5])
end

% Finally label the two panels as A and B
annotation('textbox',[0.00 0.98 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.00 0.51 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
