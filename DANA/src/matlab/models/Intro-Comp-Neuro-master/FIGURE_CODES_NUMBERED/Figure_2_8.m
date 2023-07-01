% Figure_2_8.m
% Adaptive Exponential Leaky Integrate and Fire Model
% a 2-variable model than can reproduce many types of neural behavior
% see Naud, Marcille, Clopath, Gerstner, Biol. Cybern. 2006
% this code loops through different applied currents
%
% This code was used to produce Figure 2.8 in the textbook:
% "An Introductory Course in Computational Neuroscience" 
% by Paul Miller, Brandeis University 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;                     % Begin codes by clearing memory

%% List of Cell Parameters

G_L = 10e-9;               % Leak conductance (S)
C = 100e-12;               % Capacitance (F)
E_L = -70e-3;              % Leak potential (V)
V_Thresh = -50e-3;         % Threshold potential (V)
V_Reset = -80e-3;          % Reset potential (V)
deltaT = 2e-3;             % Threshold shift factor (V)
tau_sra = 200e-3;          % Adaptation time constant (s)
a = 2e-9;                  % adaptation recovery (S)
b = 0.02e-9;               % adaptation strength (A)

I0 = 0e-9;                 % Baseline current (A)

Vmax = 50e-3;              % level of voltage to detect a spike

%% Simulation set-up

dt = 2e-6;                  % time-step in sec
tmax = 5;                   % maximum time in sec
tvector = 0:dt:tmax;        % vector of all the time points

ton = 0;                    % time to switch on current step
toff = tmax;                % time to switch off current step
non = round(ton/dt)+1;      % index of time vector to switch on
noff = round(toff/dt);      % index of time vector to switch off
I = I0*ones(size(tvector)); % baseline current added to all time points

Iappvec = [0.15:0.005:0.3]*1e-9;    % list of applied currents

initialrate = zeros(size(Iappvec)); % array to store 1/(first ISI)
finalrate = zeros(size(Iappvec));   % array to store 1/(final ISI)
singlespike = zeros(size(Iappvec)); % array to store "1" for only 1 spike
meanV = zeros(size(Iappvec));

trial = 0;                          % count which trial we are on
for Iapp = Iappvec;                 % loop through applied currents
    trial = trial+1;                % add one to the trial counter
    I(non:noff) = Iapp;             % update with new applied current
    
    
    v = zeros(size(tvector));       % initialize voltage array
    v(1) = E_L;                      % set value of initial membrane potential
    I_sra = zeros(size(tvector));   % initialize adaptation current
    spikes = zeros(size(tvector));  % initialize vector to store spikes
    
    for j = 1:length(tvector)-1     % simulation for all time points
        
        if ( v(j) > Vmax )              % if there is a spike
            v(j) = V_Reset;             % reset the voltage
            I_sra(j) = I_sra(j) + b;    % increase the adaptation current by b
            spikes(j) = 1;              % record the spike
        end
        
        % next line integrates the voltage over time, first part is like LIF
        % second part is an exponential spiking term
        % third part includes adaptation
        v(j+1) = v(j) + dt*( G_L*(E_L-v(j) + deltaT*exp((v(j)-V_Thresh)/deltaT) ) ...
            - I_sra(j) + I(j))/C;
        
        % next line decys the adaptation toward a steady state in between spikes
        I_sra(j+1) = I_sra(j) + dt*( a*(v(j)-E_L) - I_sra(j) )/tau_sra;
        
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
    
    meanV(trial) = mean(v);
end

%% Now plot the results
% First set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf;                                % Clear the figure
hold on;                            % Allow many plots on same graph
plot(1e9*Iappvec,finalrate,'k')     % Plot the final rate as a line
ISIindices = find(initialrate);     % Find points where a first ISI exists
% Then plot with an 'o' just those points where there is a first ISI
plot(1e9*Iappvec(ISIindices),initialrate(ISIindices),'ok')

ISIindices = find(singlespike);     % Find points with just one spike (no ISI)
% Then plot with a '*' on the x-axis those values that produce one spike
plot(1e9*Iappvec(ISIindices),0*singlespike(ISIindices),'*k')
xlabel('Iapp (nA)')                 % Label x-axis
ylabel('Spike Rate (Hz)')           % Label y-axis
% Now a legend labels the lines/symbols in the order of plot commands
% The legend can be moved within the figure by 'click' and 'drag'
legend('Final Rate',  '1/ISI(1)' , 'Single spike') 

