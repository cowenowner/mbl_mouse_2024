% Figure_4_12.m
% This model contains a T-type Calcium current to generate a
% post-inhibitory rebound as a model of thalamic relay cells.
% The code will step from a hyperpolarizing applied current in 5 increments
% of decreasing hyperpolarization to depolarization.
%
% This code is used to produce Figure 4.12 in the textbook 
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;                  % clear memory of prior variables

%% Simulation parameters for Fig. 4.12
dt = 0.00001;       % small time-step when action potentials are simulated
istart = 0.5;       % time applied current starts
I0=-0.1e-9;         % Initial hyperpolarizing current
Istep= 0.05e-9;     % magnitude of each step in the current
Nsteps = 6;         % Number of values of applied current
steplength = 0.25;  % Duration of each step of constant current
tmax = Nsteps*steplength;   % Total simulation time
nsteplength = round(steplength/dt); % No. of time points per current step
t=0:dt:tmax;        % time vector

%% Parameters for the Thalamic rebound model used here
E_L = -0.070;       % leak reversal potential
E_Na = 0.055;       % reversal for sodium channels
E_K = -0.090;       % reversal for potassium channels
E_Ca = 0.120        % reversal potential for Ca current

G_L = 11e-9;        % leak conductance in Siemens 
G_Na = 2.5e-6;      % sodium conductance
G_K = 1.5e-6;       % potassium conductance
G_CaT = 0.19e-6;    % T-type calcium conductance

Cm = 0.1e-9;        % membrane capacitance in Farads 

%% Initialize variables used in the simulation
I_L= zeros(size(t));    % to store leak current
I_Na= zeros(size(t));   % to store sodium current
I_K= zeros(size(t));    % to store potassium current
I_CaT = zeros(size(t)); % to store T-type calcium current

V=zeros(size(t));   % membrane potential vector
V(1) = -0.078;      % initialize membrane potential
n=zeros(size(t));   % n: potassium activation gating variable
n(1) = 0.025;       % initialize near steady state
m=zeros(size(t));   % m: sodium activation gating variable
m(1) = 0.005;       % initialize near steady state
h=zeros(size(t));   % h: sodim inactivation gating variplot(t,V)able
h(1) = 0.6;         % initialize near steady state

mca=zeros(size(t)); % CaT current activation gating variable
mca(1) = 0.025;     % initialize near steady state   
hca=zeros(size(t)); % CaT current inactivation gating variable
hca(1) = 0.6;       % initialize near steady state

Iapp=zeros(size(t)); % Applied current vector
for step = 1:Nsteps;                    % Loop through current steps
    istart = (step-1)*nsteplength+1;    % Index of current onset
    istop = step*nsteplength;           % Index of current offset
    Iapp(istart:istop) = I0+(step-1)*Istep; % Set applied current
end

Itot=zeros(size(t)); % in case we want to plot and look at the total current

%% Commence the simulation through time

for i = 2:length(t); % now see how things change through time
    Vm = V(i-1); 
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta. Written out from Dayan/Abbott, units are 1/sec.
    if ( Vm == -35 ) 
        alpha_m = 1e3;
    else 
        alpha_m = 1e5*(Vm+0.035)/(1-exp(-100*(Vm+0.035)));
    end
    beta_m = 4000*exp(-(Vm+0.060)/0.018);

    % Now sodium inactivation rate constants
    alpha_h = 350*exp(-50*(Vm+0.058));
    beta_h = 5000/(1+exp(-100*(Vm+0.028)));
    
    % Now potassium activation rate constants (the "if" prevents a divide
    % by zero)
    if ( Vm == -0.034 ) 
       alpha_n = 500;
    else
        alpha_n = 5e4*(Vm+0.034)/(1-exp(-100*(Vm+0.034)));
    end
    beta_n = 625*exp(-12.5*(Vm+0.044));
     
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.   
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);      % time constant converted from ms to sec
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);      % time constant converted from ms to sec
    n_inf = alpha_n/(alpha_n+beta_n);   
    
    % for the Ca_T current gating variables are given by formulae for the 
    % steady states and time constants:    
    mca_inf = 1/(1+exp(-(Vm+0.052)/0.0074));    % Ca_T activation
    hca_inf = 1/(1+exp(500*(Vm+0.076)));        % Ca_T inactivation
    if ( Vm < -80 ) 
        tau_hca = 1e-3*exp(15*(Vm+0.467));
    else
        tau_hca = 1e-3*(28+exp(-(Vm+0.022)/0.0105));
    end

    m(i) = m_inf;    % Update m, assuming time constant is neglible.
    
    h(i) = h_inf - (h_inf-h(i-1))*exp(-dt/tau_h);    % Update h
    
    n(i) = n_inf - (n_inf-n(i-1))*exp(-dt/tau_n);    % Update n
        
    mca(i) = mca_inf;                           % Update mca instantaneously
    hca(i) = hca_inf - (hca_inf-hca(i-1))*exp(-dt/tau_hca); % update hca
    
    G_Na_now = G_Na*m(i)*m(i)*m(i)*h(i);    % sodium conductance
    I_Na(i-1) = G_Na_now*(E_Na-V(i-1));     % sodium current
    
    G_K_now = G_K*n(i)*n(i)*n(i)*n(i);      % potassium conductance
    I_K(i-1) = G_K_now*(E_K-V(i-1));        % potassium current
    
    G_CaT_now = G_CaT*mca(i)*mca(i)*hca(i); % T-type calcium conductance
    I_CaT(i-1) = G_CaT_now*(E_Ca-V(i-1));   % Calcium T-type current
        
    I_L(i-1) = G_L*(E_L-V(i-1));            % Leak current

    Itot(i-1) = I_L(i-1)+I_Na(i-1)+I_K(i-1) ...
                +I_CaT(i-1) +Iapp(i-1); % total current is sum of leak + active channels + applied current
     
    G_Tot = G_L+G_Na_now+G_K_now+G_CaT_now; % Total conductance

    % V_inf is steady state voltage given all conductances and reversals
    V_inf = (G_L*E_L + G_Na_now*E_Na + G_K_now*E_K  + G_CaT_now*E_Ca+Iapp(i-1))/G_Tot;
        
    % Membrane potential update is via the exponential Euler method
    V(i) = V_inf - (V_inf-V(i-1))*exp(-dt*G_Tot/Cm);  

    
end

%% Now plot the results for Figure 3.12
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf
% First plot the applied current vs time
subplot('Position',[0.14 0.735 0.82 0.24])
plot(t,Iapp*1e9,'k')
ylabel('I_{app} (nA)')
axis([0 tmax -0.12 0.17])
set(gca,'XTick',[0 0.25 0.5 0.75 1 1.25 1.5])

% Next plot the membrane potential vs time
subplot('Position',[0.14 0.425 0.82 0.24])
plot(t,V*1000,'k')
ylabel('V_{m} (mV)')
axis([0 tmax -85 50])
set(gca,'XTick',[0 0.25 0.5 0.75 1 1.25 1.5])

% Finally plot the T-type calcium inactivation variable
subplot('Position',[0.14 0.115 0.82 0.24])
plot(t,hca,'k')
xlabel('Time (sec)')
ylabel('h_{CaT} ')
axis([0 tmax 0 1])
set(gca,'XTick',[0 0.25 0.5 0.75 1 1.25 1.5])

annotation('textbox',[0 0.97 0.03 0.04],'String','A','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.66 0.03 0.04],'String','B','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.35 0.03 0.04],'String','C','LineStyle','none','FontSize',16,'FontWeight','Bold')


    