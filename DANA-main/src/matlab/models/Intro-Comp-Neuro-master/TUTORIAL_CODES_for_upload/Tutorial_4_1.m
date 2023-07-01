% Tutorial_4_1.m 
%
% Hodgkin-Huxley model (new, modern units and with absolute rather than 
% specific parameters assuming a mambrane surface area of 0.1 mm^2)
%
% This code provides a solution to Tutorial 4.1 in the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University
% (February 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

part = 'e';         % Options are 'b', 'c', 'd', 'e', 'f' for question part

dt = 2e-8;          % Time step (s) must be extremely small
tmax=0.35;          % Maximum simulation time
t=0:dt:tmax;        % Time vector

%% Set up parameters and initial conditions

V_L = -0.060;       % Leak reversal potential (V)
E_Na = 0.045;       % Reversal for sodium channels (V)
E_K = -0.082;       % Reversal for potassium channels (V)
V0 = -0.065;

G_L = 30e-9;        % Leak conductance (S)
G_Na = 12e-6;       % Sodium conductance (S)
G_K = 3.6e-6;       % Potassium conductance (S)

Cm = 100e-12;       % Membrane capacitance (F)

istart = 100e-3;    % Default time applied current starts
ilength= 5e-3;      % Default length of applied current pulse
Ibase = 0e-9;       % Default baseline current before/after pulse
Npulses = 1;        % Default number of current pulses
pulsesep = 20e-3;   % Default separation between current pulses
V0 = -0.065;        % Default initial condition for V
m0 = 0.05;          % Default initial condition for m
h0 = 0.5;           % Default initial condition for h
n0 = 0.35;          % Default initial condition for n

switch part;        % Change some parameters or initial conditions with question part
    case 'b'        
        ilength= 100e-3;    % length of applied current pulse in Part b
        Ie = 0.22e-9;       % applied current during pulse in part b    
    case 'c'
        Npulses = 10;       % 10 pulses in part c
        Ie = 0.22e-9;       % applied current for each pulse in part c
        pulsesep = 18e-3;   % time between pulses in part c
    case 'd'
        Npulses = 10;       % 10 pulses in part d
        Ibase = 0.6e-9;     % large baseline current in part d
        Ie = 0e-9;          % pulse current is below baseline in part d
    case 'e'                    
        Ibase = 0.65e-9;    % baseline current in part e
        Ie = 1e-9;          % pulsed current in part e
    case 'f'        
        Ibase = 0.65e-9;    % baseline current in part f (same as e)
        Ie = 1e-9;          % pulsed current in part f (same as e)
        m0 = 0;             % initial condition for m in part f
        h0 = 0;             % initial condition for h in part f
        n0 = 0;             % initial condition for n in part f
    otherwise
        fprintf('Variable part must be b-f')
end

%% Set up the applied current vector
Iapp=Ibase*ones(size(t));       % Initialize current vector at baseline 

for pulse = 1:Npulses;          % For each current pulse
    pulsestart = istart + (pulse-1)*pulsesep;   % Onset time of pulse
    pulsestop = pulsestart + ilength;           % Offset time of pulse
    
    % make applied current a new value for duration of current pulse
    for i=round(pulsestart/dt)+1:round(pulsestop/dt);
        Iapp(i) = Ie;
    end
end

%% Now initialize all variables in the simulation

V=zeros(size(t));       % voltage vector
V(1) = V0;              % set the inititial value of voltage

n=zeros(size(t));       % n: potassium activation gating variable
n(1) = n0;               % initialize as zero
m=zeros(size(t));       % m: sodium activation gating variable
m(1) = m0;               % initialize as zero
h=zeros(size(t));       % h: sodim inactivation gating variable
h(1) = h0;               % initialize as zero

Itot=zeros(size(t));    % to plot and look at the total current
I_Na=zeros(size(t));    % to plot and look at sodium current
I_K=zeros(size(t));     % to plot and look at potassium current
I_L=zeros(size(t));     % to plot and look at leak current

for i = 2:length(t); % now see how things change through time
    
    Vm = V(i-1);          % membrane potential for calculations
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else;                   % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.070)/0.08);     % potassium deactivation rate
        
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    
    m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
    
    h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
    
    n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
    
    I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i-1)); % total sodium current
    
    I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i-1)); % total potassium current
    
    I_L(i) = G_L*(V_L-V(i-1));    % Leak current is straightforward
    
    Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i); % total current is sum of leak + active channels + applied current
    
    V(i) = V(i-1) + Itot(i)*dt/Cm;        % Update the membrane potential, V.
    
end

%% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

hold on;

figure(1)
clf

subplot('Position',[0.15 0.6 0.8 0.36])
plot(t(10:10:end),1e9*Iapp(10:10:end),'k')
ylabel('I_{app} (nA)')
axis([0 tmax -0.5 1.05])
hold on

titlestring = ['4.1 Part  ' part]
title(titlestring);

subplot('Position',[0.15 0.1 0.8 0.36])
plot(t(10:10:end),1e3*V(10:10:end),'k');


xlabel('Time, sec')
ylabel('V_{m} (mV)')
if ( max(V) > 0 )
    axis([0 tmax -85 45])
else
    axis([0 tmax -80 -55])
end
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.5 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')
