% Tutorial_3_3.m
% Tutorial_3_3.m produces a Receiver Operating Characteristic curve for
% a noisy AELIF neuron receiving two values of input current corresponding
% to a stimulus or its absence.
% The input current has a stimulus-dependent mean and standard deviation
% that can change in each part of the tutorial. Stimulus time is reduced in
% Part 2.
% AELIF is the Adaptive Exponential Leaky Integrate and Fire Model,
% a 2-variable model than can reproduce many types of neural behavior
% see Naud, Marcille, Clopath, Gerstner, Biol. Cybern. 2006
%
% This code is a solution to Tutorial 3.3 in Chapter 3 of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;                      % Clear all variables
Ntrials = 1000;             % No. of trials to obtain probability distribution
%% List of Cell Parameters

G_L = 10e-9;               % Leak conductance (S)
C = 100e-12;               % Capacitance (F)
E_L = -70e-3;              % Leak potential (V)
V_Thresh = -50e-3;         % Threshold potential (V)
V_Reset = -80e-3;          % Reset potential (V)
deltaT = 2e-3;             % Threshold shift factor (V)
tauw = 150e-3;             % Adaptation time constant (s)
a = 2e-9;                  % adaptation recovery (S)
b = 0e-9;                  % adaptation strength (A)

Vmax = 50e-3;              % level of voltage to detect a spike

I0 = zeros(1,2);
sigma_I = zeros(1,2);
for question = 1:3
    %% Simulation set-up
    
    dt = 5e-6;                          % dt in sec
    switch question
        case 1;
            tmax = 0.5;                 % Duration of each trial
            I0(1) = 0.0e-9;             % Mean applied current for no stimulus
            I0(2) = 0.1e-9;             % Mean applied current with stimulus present
            sigma_I(1)= 0.02e-9;        % S.D. of input current with no stimulus
            sigma_I(2)= 0.02e-9;        % S.D. of input current with a stimulus
        case 2;
            tmax = 0.2;                 % Duration of each trial
            I0(1) = 0.0e-9;             % Mean applied current for no stimulus
            I0(2) = 0.1e-9;             % Mean applied current with stimulus present
            sigma_I(1)= 0.02e-9;        % S.D. of input current with no stimulus
            sigma_I(2)= 0.02e-9;        % S.D. of input current with a stimulus
        case 3;
            tmax = 0.5;                 % Duration of each trial
            I0(1) = 0.0e-9;             % Mean applied current for no stimulus
            I0(2) = 0.3e-9;             % Mean applied current with stimulus present
            sigma_I(1)= 0.05e-9;        % S.D. of input current with no stimulus
            sigma_I(2)= 0.005e-9;       % S.D. of input current with a stimulus
    end
    
    tvector = 0:dt:tmax;        % Vector of time points
    
    % Nspikes stores number of spikes of each trial in columns, with the two
    % rows corresponding to the two values of the stimulus (on/off)
    Nspikes = zeros(2,Ntrials);
    
    for stimulus = 1:2;         % Stimulus can be on or off
        for trial = 1:Ntrials;  % Loop through all the trials for each stimulus
            
            % Iapp contains a stimulus-dependent constant component plus a
            % noise component that could also be stimulus-dependent
            Iapp = I0(stimulus)*ones(size(tvector)) ...
                + sigma_I(stimulus)*randn(size(tvector))/sqrt(dt);
            
            v = zeros(size(tvector));       % initialize voltage for the trial
            v(1) = E_L;                     % set value at t=0
            w = zeros(size(tvector));       % initialize adaptation variable
            spikes = zeros(size(tvector));  % initially no spikes for the trial
            
            for j = 1:length(tvector)-1     % simulation up to tmax
                
                if ( v(j) > Vmax )          % if there is a spike
                    v(j) = V_Reset;         % reset the voltage
                    w(j) = w(j) + b;        % increase the adaptation variable by b
                    spikes(j) = 1;          % record the time-bin (j) of the spike
                end
                
                % next line integrates the voltage over time, first part is like LIF
                % second part is an exponential spiking term
                % third part includes adaptation
                v(j+1) = v(j) + dt*( G_L*(E_L-v(j) + deltaT*exp((v(j)-V_Thresh)/deltaT) ) ...
                    - w(j) + Iapp(j))/C;
                
                % next line decays the adaptation toward a steady state in between spikes
                w(j+1) = w(j) + dt*( a*(v(j)-E_L) - w(j) )/tauw;
                
            end
            Nspikes(stimulus,trial) = sum(spikes);  % Record total no. of spikes in trial
        end
    end
    
    Nmax = max(max(Nspikes));                       % Max number of spikes in any trial
    counts1 = histcounts(Nspikes(1,:),[0:Nmax]);          % Histogram of spike counts (range 0 to Nmax)
    counts2 = histcounts(Nspikes(2,:),[0:Nmax]);          % Histogram of spike counts (range 0 to Nmax)
    
    % The cumulative sum of the histogram is taken as "cumsum", giving the
    % number of trials with that number or fewer. Division by Ntrials converts
    % the number of trials to the fraction of trials with that number or fewer.
    % Subtraction of this from 1 creates the fraction greater than each value.
    frac_counts1 = cumsum(counts1,'reverse')/Ntrials;         % Fraction of counts greater than each value
    frac_counts2 = cumsum(counts2,'reverse')/Ntrials;         % Fraction of counts greater than each value
    
    %% Set up the plotting parameters
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(question)        % New figure for each of questions 1-3
    clf
    subplot('Position',[0.25 0.74 0.7 0.23])        % Top of three figures
    stairs([0.5:Nmax],counts1,'k','LineWidth',3);     % Histogram with no stimulus
    hold on
    stairs([0.5:Nmax],counts2,'r','LineWidth',3);     % Histogram with stimulus
    axis([0 Nmax+1 0 max([counts1,counts2])*1.05])
    xlabel('Spike count')
    ylabel('No. of trials')
    
    subplot('Position',[0.25 0.41 0.7 0.23])        % Second of three figures
    stairs([0.5:Nmax],frac_counts1,'k','LineWidth',3) % 1 - cumulative probability
    hold on
    stairs([0.5:Nmax],frac_counts2,'r','LineWidth',3) % 1 - cumulative probability
    axis([0 Nmax+1 0 1])
    xlabel('Spike count')
    ylabel('Prob. greater')
    
    subplot('Position',[0.25 0.08 0.7 0.23])        % Third of three figures
    plot(frac_counts1,frac_counts2,'k')                 % ROC curve
    xlabel('Prob. of False +ve')
    ylabel('Prob. of True +ve')
    
    
    annotation('textbox',[0 0.98 0.03 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
    annotation('textbox',[0 0.65 0.03 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
    annotation('textbox',[0 0.32 0.03 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
    drawnow
end