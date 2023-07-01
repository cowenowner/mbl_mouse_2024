% HHab_plotnewV.m 
% Plots currents and gating variables as a function of V
% This version is in Hodgkin-Huxley's original formalism, 
% where the reversal potential is defined to be 0 and what we 
% now call depolarization (via an inward current) corresponds to 
% a negative potential difference.

% Note: The plots have been switched to plot out in "modern units" 
% of voltage by replacing "V" with "-70-V" (in mV). 
%
% This code is used to produce Figures 4.3 and 4.4 in the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;              % clear all prior variables and parameters from memory
dV = 0.1;           % step size of membrane potential 
Vmin = -100;        % minimum of range for the calculation
Vmax = 50;          % maximum of range for the calculation
V = Vmin:dV:Vmax;   % vector of voltage values

V_L = -10.613;      % leak reversal potential (mV in old units)
E_Na = -115;        % reversal for sodium channels (mV in old units)
E_K = 12;           % reversal for potassium channels (mV in old units)

g_L = 0.3;          % specific leak conductance (M-ohms /cm^2)
g_Na = 120;         % specific sodium conductance (M-ohms /cm^2)
g_K = 36;           % specific potassium conductance (M-ohms /cm^2)

cm = 1;             % specific membrane capacitance (nF /cm^2)
tau = cm/g_L;       % membrane time constant (ms)

n=zeros(size(V));   % n: potassium activation gating variable
m=zeros(size(V));   % m: sodium activation gating variable
h=zeros(size(V));   % h: sodim inactivation gating variable
% The set of alpha and beta are the rate constants for each of the three
% gating variables
alpha_n = zeros(size(V));   % potassium activation opening rate
beta_n = zeros(size(V));    % potassium activation closing rate
alpha_m = zeros(size(V));   % sodium activation opening rate
beta_m = zeros(size(V));    % sodium activation closing rate
alpha_h = zeros(size(V));   % sodium inactivation opening rate
beta_h = zeros(size(V));    % sodium inactivation closing rate

I_L=zeros(size(V));         % vector to store the leak current
I_Na=zeros(size(V));        % vector to store the sodium current
I_K=zeros(size(V));         % vector to store the potassium current
Itot=zeros(size(V));        % vector to store the total current

for i = 1:length(V); % now calculate steady states as a function of V

    I_L(i) = g_L*(V_L-V(i));    % leak current

    Vm = V(i);          % allows for easy switch of units

    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta. Written out from Dayan/Abbott, units are 1/ms.
    if ( Vm == -25 )
        alpha_m(i) = 0.1/0.1;
    else
        alpha_m(i) = (0.1*(Vm+25))/(exp(0.1*(Vm+25))-1);
    end
    beta_m(i) = 4*exp(Vm/18);
    alpha_h(i) = 0.07*exp(Vm/20);
    beta_h(i) = 1/(1+exp((Vm+30)/10));
    if ( Vm == -10)
        alpha_n(i) = 0.01/0.1;
    else
        alpha_n(i) = (0.01*(Vm+10))/(exp(0.1*(Vm+10))-1);
    end
    beta_n(i) = 0.125*exp((Vm)/80);

    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    tau_m(i) = 1/(alpha_m(i)+beta_m(i));      
    m(i) = alpha_m(i)/(alpha_m(i)+beta_m(i));

    tau_h(i) = 1/(alpha_h(i)+beta_h(i));     
    h(i) = alpha_h(i)/(alpha_h(i)+beta_h(i));

    tau_n(i) = 1/(alpha_n(i)+beta_n(i));      
    n(i) = alpha_n(i)/(alpha_n(i)+beta_n(i));

    % Now calculate currents from gating variables
    I_Na(i) = g_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i)); % total sodium current
    I_K(i) = g_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i)); % total potassium current
    Itot(i) = I_L(i)+I_Na(i)+I_K(i); % total current is sum of leak + active channels + applied current


end

%% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

% figure(1) compares the different steady-state currents as a function of V
figure(1)
clf
plot(-70-V,Itot);
hold on;
plot(-70-V,I_Na,'g');
plot(-70-V,I_K,'r');
plot(-70-V,I_L,'m');
xlabel('Voltage, mV ')
ylabel('Current, nA')
legend('Total', 'Na', 'K', 'Leak')

% figure(2) compares the steady states of the gating variables as a
% function of V
% This is Figure 4.4 of the textbook
figure(2);
clf
plot(-70-V,m,'k');
hold on;
plot(-70-V,h,'k--');
plot(-70-V,n,'k:');
xlabel('Voltage, mV')
ylabel('Gating variable steady state')
legend('m^{\infty }', 'h^{\infty }', 'n^{\infty }')
axis([-100 20 0 1])

% figure(3) compares the rate constants as a function of V (3 panels)
% This is Figure 4.3 of the textbook
figure(3);
clf
% First panel is sodium activation
subplot('Position',[0.108 0.15 0.24 0.8])
plot(-70-V,alpha_m,'k');
hold on;

plot(-70-V,beta_m,'k:');

xlabel('Voltage (mV)')
ylabel('rate constants (ms^{-1})')
legend('\alpha_m' , '\beta_m')
axis([-100 20 0 20])

% Second panel is potassium activation
subplot('Position',[0.42 0.15 0.24 0.8])
plot(-70-V,alpha_n,'k');
hold on;

plot(-70-V,beta_n,'k:');

xlabel('Voltage (mV)')
legend('\alpha_n' , '\beta_n')
axis([-100 20 0 1])

% Third panel is sodium inactivation
subplot('Position',[0.72 0.15 0.24 0.8])
plot(-70-V,alpha_h,'k');
hold on;

plot(-70-V,beta_h,'k:');
xlabel('Voltage (mV)')
legend('\alpha_h' , '\beta_h')
axis([-100 20 0 1])


