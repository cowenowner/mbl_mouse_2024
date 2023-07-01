% Figure_4_5.m
%
% This model is the Hodgkin-Huxley model simulated in old units, but
% plotted in new units (V switched to -70mV - V).
%
% This code is used to produce Figures 4.2 and 4.5 in the textbook,
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;                  % Clear all variables and parameters from memory
dt = 0.0002;            % Small time-step of 0.2 micro-sec for these simulations 
tmax=200;               % Simulate for 200ms (Time measured in ms)

istart = 50;            % time applied current starts (ms)
ilength=100;            % length of applied current pulse (ms)
Ibase = 0e-7;           % base current density (nA/cm^2) outside of pulse
Ie= 7;                  % extra current density (nA/cm^2) during pulse

V_L = -10.613;          % leak reversal potential (mV, old units)
E_Na = -115;            % reversal for sodium channels  (mV, old units)
E_K = 12;               % reversal for potassium channels  (mV, old units)

g_L = 0.3;              % specific leak conductance (M-ohms /cm^2)
g_Na = 120;             % specific sodium conductance (M-ohms /cm^2)
g_K = 36;               % specific potassium conductance (M-ohms /cm^2)

cm = 1;                 % specific membrane capacitance (nF / cm^2)

t=0:dt:tmax;            % time vector
V=zeros(size(t));       % voltage vector

Iapp=Ibase*ones(size(t)); % Applied current, 
for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
    Iapp(i) = Ie;
end

V(1) = V_L;          % set the inititial value of membrane potential

n=zeros(size(t));   % n: potassium activation gating variable
n(1) = 0.35;        % start off near steady state when V is at leak potential
m=zeros(size(t));   % m: sodium activation gating variable
m(1) = 0.05;        % start off near steady state when V is at leak potential
h=zeros(size(t));   % h: sodim inactivation gating variable
h(1) = 0.75;        % start off near steady state when V is at leak potential

Itot=zeros(size(t));    % to plot and look at the total current
I_Na=zeros(size(t));    % to plot and look at the sodium current
I_K=zeros(size(t));     % to plot and look at the potassium current
I_L=zeros(size(t));     % to plot and look at the leak current

for i = 1:length(t)-1; % now see how things change through time
    I_L(i) = g_L*(V_L-V(i));        % calculate leak current
        
    Vm = V(i);         % just to allow changes of units if necessary here
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta. Written out from Dayan/Abbott textbook, units are 1/ms.  
    if ( Vm == -25 )
        alpha_m = 0.1/0.1;
    else
        alpha_m = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);
    end
    beta_m = 4*exp(Vm/18);
    alpha_h = 0.07*exp(Vm/20);
    beta_h = 1/(1+exp((Vm+30)/10));
    if ( Vm == -10)
        alpha_n = 0.01/0.1;
    else
        alpha_n = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);
    end
    beta_n = 0.125*exp((Vm)/80);

    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.

    tau_m = 1/(alpha_m+beta_m);      
    m_inf = alpha_m/(alpha_m+beta_m);

    tau_h = 1/(alpha_h+beta_h);     
    h_inf = alpha_h/(alpha_h+beta_h);

    tau_n = 1/(alpha_n+beta_n);      
    n_inf = alpha_n/(alpha_n+beta_n);

    
    if ( i > 1 ) 
        m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
    
        h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
    
        n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
    end
    
    I_Na(i) = g_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % total sodium current
    
    I_K(i) = g_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i)); % total potassium current
    
    Itot(i) = I_L(i)+I_Na(i)+I_K(i)-Iapp(i); % total current is sum of leak + active channels + applied current
    
    V(i+1) = V(i) + Itot(i)*dt/cm;        % Update the membrane potential, V.
    
end

%% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

hold on;

% figure(1) plots V in mV and applied current in nA/cm^2 on same graph
figure(1)
clf
plot(t,-70-V,'k');
hold on
plot(t,Iapp,'k--')

legend('V_{m}','I_{app}')

xlabel('Time, sec')
ylabel('V_{m} (mV)  /  I_{app} (nA)')

% figure(2) plots magnitude of sodium and potassium currents on top of the
% spike waveform
figure(2)
clf
plot(t,(E_K-V)*5);
hold on
plot(t,-I_Na,'g')
plot(t,I_K,'r')
axis([66 76 0 800])

legend('scaled V', 'I_{Na}', '-I_K')

xlabel('time, sec')

%% figure(3) becomes Figure 4.5 in the textbook
% figure(3) plots the three gating variables on the same curve around a
% single spike
figure(3)
clf
subplot('Position',[0.13 0.13 0.8 0.8])
plot(t,m,'k:')
hold on
plot(t,h,'k--')
plot(t,n,'k')
axis([66 76 0 1])

legend('m', 'h', 'n')
xlabel('Time (ms)')
ylabel('Gating variable')

%% figure(4) becomes Figure 4.2 in the textbook
figure(4)
clf
plot(t,-70-V,'k');
hold on
plot(t,(-70-E_K)*ones(size(t)),':k')
plot(t,(-70-E_Na)*ones(size(t)),':k')
axis([66 76 -90 50])

ylabel('Membrane Potential (mV)')
xlabel('Time (ms)')

