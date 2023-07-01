% LIF_model.m
% Simulates a leaky integrate-and-fire neuron, which receives an applied
% current pulse of three different amplitudes.
% This code was used to produce Figure 2.3 in the textbook:
% "An Introductory Course in Computational Neuroscience"
% by Paul Miller, Brandeis University 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the default plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf;

%% Simulation parameters
dt = 0.0001;                % time-step
t = 0:dt:0.5;               % vector of time-points
ton = 0.15;                 % time to begin applied current (onset)
toff = 0.35;                % time to end applied current (offset)
non = round(ton/dt);        % time-point index of current onset
noff = round(toff/dt);      % time-point index of current offset

%% Parameters for the LIF neuron
tau = 0.010;                % membrane time constant
E_L = -0.070;               % leak potential (also resting potential)
Vth = -0.050;               % threshold potential (to produce spike)
Vreset = -0.080;            % reset potential (post-spike)
Cm = 100e-12;               % total membrane capacitance
G_L = Cm/tau;               % total membrane conductance (leak conductance)

%% Now simulate trials, each with a different applied current
%
Iapp = [180e-12 210e-12 240e-12];   % values of applied current steps
Ntrials = length(Iapp);             % number of different trials
for trial = 1:Ntrials;              % loop through different trials
    I = zeros(size(t));             % vector for current at each time-point
    I(non:noff) = Iapp(trial);      % add the applied current for the trial
    V = E_L*ones(size(t));          % initialize the membrane potential vector
    spikes = zeros(size(t));        % initialize a vector to record spikes
    
    for i = 2:length(t);            % loop through all time points
        % next line: Forward Euler method to update membrane potential
        % see Eq. 1.9 in the book
        V(i) = V(i-1) + dt*(I(i) +G_L*(E_L-V(i-1)))/Cm;
        if V(i) > Vth;              % if potential is above threshold
            spikes(i) = 1;          % record the spike at that time-point
            V(i) = Vreset;          % reset the potential
        end;
    end;                            % end the loop & go to next time-point
    
    % Now generate a series of subplots to position the subfigures
    subplot('position',[0.31*trial-0.18 0.74 0.22 0.22])
    
    plot(t,I*1e9,'k');              % Plot current (in nA) against time
    if ( trial == 1 )
        ylabel('I_{app} (nA)')
    end
    % The following line adds a title, which uses the command "strcat" to
    % combine the value of current in that trial with "nA" to make a
    % complete title.
    % The command "num2str" converts the number representing the applied
    % current into a string, which is something that can be printed.
    title(strcat(num2str(1e9*Iapp(trial)),'nA'))
    
    set(gca,'XTick',[0 0.25 0.5])
    
    axis([0 0.5 0 0.3])         % Sets x-axis then y-axis ranges
    
    % Next subfigure
    subplot('position',[0.31*trial-0.18 0.43 0.22 0.22])
    plot(t,1000*V,'k');                      % Plot membrane potential vs time
    if ( trial == 1 )
        ylabel('V_m (mV)')               % Label y-axis
    end
    set(gca,'XTick',[0 0.25 0.5])
    axis([0 0.5 1000*(Vreset-0.005) 1000*(Vth+0.005)])    % Set the axes
    
    subplot('position',[0.31*trial-0.18 0.12 0.22 0.22])
    plot(t,spikes,'k')          % Plots a line from 0 to 1 at each spike
    xlabel('Time (sec)')        % Label the x-axis
    if ( trial == 1)
        ylabel('Spikes')        % Label the y-axis
    end
    set(gca,'XTick',[0 0.25 0.5])
    set(gca,'YTick',[0 1])
    axis([0 0.5 0 1])           % Set the ranges of x- and y-axes
    
end;                            % Loop to the next trial with new current

annotation('textbox',[0.00 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.37 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.68 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C')

