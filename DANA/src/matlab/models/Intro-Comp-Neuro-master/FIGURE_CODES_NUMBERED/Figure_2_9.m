% Figure_2_9.m
% Variant of the Adaptive Exponential Leaky Integrate and Fire Model,
% a 2-variable model than can reproduce many types of neural behavior
% see Naud, Marcille, Clopath, Gerstner, Biol. Cybern. 2006
% This code loops through different applied currents.
% The adaptation term is altered in three ways:
% 1) It is a conductance rather than a current.
% 2) It has a time constant of 0.5ms, so acts as a refractory,
% rate-limiting term rather than a standard spike-rate adaptation term.
% 3) A tan function of voltage is used, so the conductance increases to
% very high levels if the membrane potential is high. The tan function is
% so steep that it can counteract excessively strong input currents that
% would otherwise cause the mean membrane potential to constantly grow with
% input. This prevents the firing rates from increasing without end.
%
% The threshold term is altered as in Method 2 of Tutorial 2.2, reaching a
% peak then decaying back to baseline following each spike.
% The reset potential is set to the baseline of the threshold, so plays a 
% greater role at higher rates when threshold is reached before it has 
% returned to baseline. 
% A dip in membrane potential following reset arises from the refractory
% conductance.
%
% This code was used to produce Figure 2.9 in the text book:
% An Introductory Course in Computational Neuroscience 
% by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;                     % Begin codes by clearing memory

%% List of Cell Parameters

G_L = 10e-9;                % Leak conductance (S)
C = 100e-12;                % Capacitance (F)
E_L = -75e-3;               % Leak potential (V)
E_K = -80e-3;               % Leak potential (V)
V_Thresh_base = -60e-3;     % Baseline of threshold potential (V)
Vmax = 150e-3;              % level of voltage to detect a spike
V_Thresh_peak = 4*Vmax;     % Maximum of threshold, post-spike
deltaT = 10e-3;             % Threshold shift factor (V)
V_Reset = -60e-3;           % Reset potential (V)

tau_g = 0.5e-3;          % Adaptation time constant (s)
tau_vth = 2e-3;

%% Setting a and b to zero converts AELIF to ELIF model

a = 10e-9;                      % adaptation recovery (S)
b = 60e-9;                      % adaptation strength (A)

I0 = 0e-9;                      % Baseline current (A)


%% Simulation set-up

dt = 5e-7;                % time-step in sec
tmax = 1;                 % maximum time in sec
tvector = 0:dt:tmax;        % vector of all the time points

ton = 0;                    % time to switch on current step
toff = tmax;                % time to switch off current step
non = round(ton/dt)+1;      % index of time vector to switch on
noff = round(toff/dt);      % index of time vector to switch off
I = I0*ones(size(tvector)); % baseline current added to all time points

Iappvec = [0:10:1000]*1e-12;    % list of applied currents

initialrate = zeros(size(Iappvec)); % array to store 1/(first ISI)
finalrate = zeros(size(Iappvec));   % array to store 1/(final ISI)
singlespike = zeros(size(Iappvec)); % array to store "1" for only 1 spike
meanV = zeros(size(Iappvec));

%% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf;                                % Clear the figure

%% Loop through trials with different applied currents

trial = 0;                          % count which trial we are on
for Iapp = Iappvec;                 % loop through applied currents
    trial = trial+1;                % add one to the trial counter
    I(non:noff) = Iapp;             % update with new applied current    
    
    v = zeros(size(tvector));       % initialize voltage array
    v(1) = E_L;                      % set value of initial membrane potential
    G_sra = zeros(size(tvector));   % initialize adaptation current
    spikes = zeros(size(tvector));  % initialize vector to store spikes
    
    inspike = 0;
    V_Thresh = V_Thresh_base*ones(size(tvector));
    
    for j = 1:length(tvector)-1     % simulation for all time points
        
        if ( v(j) > Vmax+V_Thresh(j) )      % if there is a spike
            v(j) = V_Reset;                 % reset the voltage
            G_sra(j) = G_sra(j) + b;        % increase the adaptation term 
            spikes(j) = 1;                  % record the spike
            V_Thresh(j) = V_Thresh_peak;    % increase the threshold to max
        end
        
        % next line integrates the voltage over time, first part is like LIF
        % second part is an exponential spiking term
        % third part includes adaptation 
        v(j+1) = v(j) + dt*( G_L*(E_L-v(j) + deltaT*exp((v(j)-V_Thresh(j))/deltaT) ) ...
            + (E_K-v(j))*G_sra(j) + I(j))/C;
        v(j+1) = max(v(j+1),E_K);
        
        % Next line decays the voltage threshold to baseline between spikes
        V_Thresh(j+1) = V_Thresh(j)+(V_Thresh_base-V_Thresh(j))*dt/tau_vth;
        
        % next line decays the adaptation toward a steady state that is a
        % tan function of the membrane potential (and remains at a max
        % level if the membrane potential is every above the first peak of
        % the tan function)
        Gss = a*tan(0.25*pi*min((v(j)-E_K)/(V_Thresh_base+Vmax-E_K),1.99999)); 
        G_sra(j+1) = G_sra(j) + dt*(Gss-G_sra(j))/tau_g;
    end
    
    spiketimes = dt*find(spikes);           % extract the spike times
    
    if ( length(spiketimes) > 1 )           % if there is more than 1 spike
        ISIs = diff(spiketimes);            % ISI = interval between spikes
        initialrate(trial) = 1/ISIs(1);     % inverse of first ISI
        if ( length(ISIs) > 1 )             % if there are further ISIs
            finalrate(trial) = 1/ISIs(end); % inverse of final ISI
        end
        
    else
        if ( length(spiketimes) == 1 )      % if there is only one spike
            singlespike(trial) = 1;         % record "1" for this trial
        end
    end
    
    meanV(trial) = mean(v);                 % record mean membrane potential

    if ( Iapp == 160e-12 )                  % For one applied current
        figure(1)
        subplot('Position',[0.16 0.62 0.33 0.33])                      % Plot V versus t for 200ms.
        plot(tvector(1:ceil(0.2/dt)),1e3*v(1:ceil(0.2/dt)),'k')
        xlabel('Time (sec)')
        ylabel('V_{m} (mV)')
    end
end

%% Plot the summary of results

subplot('Position',[0.65 0.62 0.33 0.33])                      % Plot V versus t for 200ms.
hold on;                            % Allow many plots on same graph
plot(1e9*Iappvec,finalrate,'k')     % Plot the final rate as a line
axis([0 1 0 120])
xlabel('I_{app} (nA)')              % Label x-axis
ylabel('Spike Rate (Hz)')           % Label y-axis

subplot('Position',[0.65 0.14 0.33 0.33])                      % Plot V versus t for 200ms.
plot(1e9*Iappvec,1e3*meanV,'k')
axis([0 1 -65 -35])
xlabel('I_{app} (nA)')

subplot('Position',[0.16 0.14 0.33 0.33])                      % Plot V versus t for 200ms.
plot(finalrate,1e3*meanV,'k')
axis([0 120 -65 -35])
xlabel('Spike rate (Hz)')
ylabel('mean of V_{m} (mV)')

annotation('textbox',[0.0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.95 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.0 0.5 0.05 0.05],'String','C','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.5 0.05 0.05],'String','D','LineStyle','none','FontSize',16,'FontWeight','Bold')
