% Figure_2_7.m
% Adaptive Exponential Leaky Integrate and Fire Model
% a 2-variable model than can reproduce many types of neural behavior
% see Naud, Marcille, Clopath, Gerstner, Biol. Cybern. 2006
%
% This code was used to produce Figure 2.7 in the textbook:
% "An Introductory Course in Computational Neuroscience" 
% by Paul Miller, Brandeis University 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% Set up the default plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
    figure(1)
    clf;


%% List of Cell Parameters

G_L = 10e-9;                % Leak conductance (S)
C = 100e-12;                % Capacitance (F) 
E_L = -70e-3;               % Leak potential (V)
V_Thresh = -50e-3;          % Threshold potential (V)
V_Reset = -80e-3;           % Reset potential (V)
deltaT = 2e-3;              % Threshold shift factor (V)
tauw = 200e-3;              % Adaptation time constant (s)
a = 2e-9;                   % Adaptation recovery (S)
b = 0.02e-9;                % Adaptation strength (A)

I0 = 0e-9;                  % Baseline current
Iapp = 0.221e-9             % Applied current step

Vmax = 50e-3;               % Level of voltage to detect and crop spikes

%% Simulation set-up

dt = 1e-6;                  % dt in sec
tmax = 3;                   % maximum time in sec
tvector = 0:dt:tmax;        % vector of time points

ton = 0.5;                  % time to add step current
toff = 2.5;                 % time to remove step current
non = round(ton/dt);        % time-point corresponding to ton
noff = round(toff/dt);      % time-point corresponding to toff
I = I0*ones(size(tvector)); % applied current initialized to baseline
I(non:noff) = Iapp;         % add the step to the applied current vector    

v = zeros(size(tvector));       % initialize membrane potential at all time-points
v(1) = E_L;                     % set initial value to be the leak potential
w = zeros(size(tvector));       % initialize adaptation variable at all time-points
spikes = zeros(size(tvector));  % intialize a vector to record spiketimes

for j = 1:length(tvector)-1     % simulation through all time-points

    if ( v(j) > Vmax )          % if there is a spike
        v(j) = V_Reset;         % reset the voltage
        w(j) = w(j) + b;        % increase the adaptation variable by b
        spikes(j) = 1;          % record the spike
    end
    
    % next line integrates the membrane potential over time, using the 
    % Forward  Euler method.
    % first term within parentheses is like the LIF model
    % second term is an exponential spiking term
    % third term includes adaptation
    v(j+1) = v(j) + dt*( G_L*(E_L-v(j) + deltaT*exp((v(j)-V_Thresh)/deltaT) ) ...
       - w(j) + I(j))/C;

   % next line decys the adaptation toward a steady state in between spikes
    w(j+1) = w(j) + dt*( a*(v(j)-E_L) - w(j) )/tauw;
    
end
            
%% Plot the results
figure(1)
subplot('Position',[0.2 0.72 0.78,0.22])       % top subpanel
plot(tvector,1e9*I,'k')     % plot input current as a function of time
ylabel('I_{app} (nA)')
axis([0 tmax 0 1.25*max(1e9*I)])

subplot('Position',[0.2 0.41 0.78,0.22])        % middle subpanel
plot(tvector,1000*v,'k')     % plot membrane potential as a function of time
axis([0 tmax -95 35])
ylabel('V_{m} (mV)') 
set(gca,'YTick',[-50 0])

subplot('Position',[0.2 0.10 0.78,0.22])       % lower subpanel
plot(tvector,1e9*w,'k')     % plot adaptation variable as a function of time
xlabel('Time (sec)')
ylabel('I_{SRA} (nA)')

annotation('textbox',[0.00 0.97 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.00 0.66 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.00 0.35 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C')
