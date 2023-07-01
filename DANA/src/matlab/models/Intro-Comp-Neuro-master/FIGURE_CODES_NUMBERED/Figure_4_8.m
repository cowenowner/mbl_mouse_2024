% Figure_4_8.m
% Simulates membrane potential as a function of time, when an outward
% current is applied to the model neuron, then released.
%
% This version is in Hodgkin-Huxley's original formalism, 
% where the reversal potential is defined to be 0 and what we 
% now call depolarization (via an inward current) corresponds to 
% a negative potential difference.

% Note: The plots have been switched to plot out in "modern units" 
% of voltage by replacing "V" with "-70-V" (in mV). 
%
% This code is used to produce Figure 4.8 in the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;                  % Clear all variables and parameters from memory
dt = 0.0002;            % Small time-step of 0.2 micro-sec for these simulations 
tmax=200;               % Simulate for 200ms (Time measured in ms)

istart = 50;            % time applied current starts (ms)
ilength=100;            % length of applied current pulse (ms)
Ibase = 0e-7;           % base current density (nA/cm^2) outside of pulse
Ie= -6;                 % extra current density (nA/cm^2) during pulse

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

V(1) = 0;           % set the initial value of membrane potential (near steady state)
n=zeros(size(t));   % n: potassium activation gating variable
n(1) = 0.35;        % start off near steady state 
m=zeros(size(t));   % m: sodium activation gating variable
m(1) = 0.05;        % start off near steady state 
h=zeros(size(t));   % h: sodim inactivation gating variable
h(1) = 0.6;         % start off near steady state 

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

% figure(1) is Figure 4.8 in the textbook
figure(1)
clf
% First panel is current versus time
subplot('Position',[0.15 0.7 0.8 0.26])
plot(t(10:10:end),Iapp(10:10:end),'k')
ylabel('I_{app} (nA)')
axis([0 tmax -6.5 0.5])
hold on

% Second panel is membrane potential versus time (shifted to modern units)
subplot('Position',[0.15 0.4 0.8 0.26])
plot(t(10:10:end),-70-V(10:10:end),'k');
axis([0 tmax -85 45])
ylabel('V_{m} (mV)')

% Third panel is gating variables versus time
subplot('Position',[0.15 0.1 0.8 0.26])
plot(t(10:10:end),m(10:10:end),'k');
hold on
plot(t(10:10:end),h(10:10:end),'k--');
plot(t(10:10:end),n(10:10:end),'k:');

axis([0 tmax 0 1])

xlabel('Time (ms)')
ylabel('Gating variable')
legend('m', 'h', 'n')

% Label the panels A, B, and C
annotation('textbox',[0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.65 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.35 0.05 0.05],'String','C','LineStyle','none','FontSize',16,'FontWeight','Bold')
