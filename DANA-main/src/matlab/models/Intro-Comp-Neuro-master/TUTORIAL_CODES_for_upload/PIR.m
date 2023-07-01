function [V, tvec] = PIR(Ibase,Istep,dt)
% PIR returns the voltage response alongside a vector of time-points for a
% simulation of duration 0.75s, split into three parts. 
% The current, Ibase, is applied to the neuron between 0 and 0.25s.
% A current of Ibase + Istep is applied between 0.25s and 0.5 s.
% A current of Ibase is applied between 0.5s and 0.75s.
% dt is the time-step used for simulation purposes (default is 2 microsec).
% 
% The model contains a T-type Calcium current to generate a
% post-inhibitory rebound as a models of thalamic relay cells.
%
% This funtion is required to complete Tutorial 4.2 of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~exist(dt) )   % If no dt was sent to the function
    dt = 2e-6;      % Use this as a default value
end
 
tmax = 0.750;       % Maximum time to simulate
tvec = 0:dt:tmax;   % Vector of time points

E_L = -0.070;       % leak reversal potential
E_Na = 0.055;       % reversal for sodium channels
E_K = -0.090;       % reversal for potassium channels
E_Ca = 0.120;        % reversal potential for Ca current

G_L = 10e-9;        % leak conductance in Siemens
G_Na = 3.6e-6;      % sodium conductance in Siemens
G_K = 1.6e-6;       % potassium conductance in Siemens
G_CaT = 0.22e-6;    % T-type calcium conductance in Siemens

Cm = 0.1e-9;        % membrane capacitance in Farads

%% Now initialize all variables
I = Ibase*ones(size(tvec));     % Baseline current
non = round(0.25/dt);           % Time-point to step up the current
noff = round(0.5/dt);           % Time-point to step the current back down
I(non:noff) = Ibase+Istep;      % Total current is baseline plus stepped amplitude

V = zeros(size(tvec));          % Array for membrane potential
V(1) = E_L;                     % Initialize at leak potential

m = zeros(size(tvec));          % Sodium activation
h = zeros(size(tvec));          % Sodium inactivation
n = zeros(size(tvec));          % Potassium activation
mca = zeros(size(tvec));        % T-type calcium activation
hca = zeros(size(tvec));        % T-type calcium inactivation

%% Now simulate through time
for i = 2:length(tvec); % now see how things change through time
    Vm = V(i-1); % For ease of repeated use in the equations
    
    % alpha_m is sodium activation rate constant
    if ( Vm == -0.035 )                            % to prevent division by zero
        alpha_m = 1e3;
    else
        alpha_m = 1e5*(Vm+0.035)./(1-exp(-100*(Vm+0.035)));
    end
    beta_m = 4000*exp(-(Vm+0.060)./0.018);      % sodium deactivation rate
    
    alpha_h = 350*exp(-50*(Vm+0.058));          % sodium inactivation rate
    beta_h = 5000./(1+exp(-100*(Vm+0.028)));    % sodium deinactivation rate
    
    if ( Vm == -0.034 )                            % prevent division by zero
        alpha_n = 500;                          % potassium activation rate constant
    else
        alpha_n = 5e4*(Vm+0.034)./(1-exp(-100*(Vm+0.034)));
    end
    beta_n = 625*exp(-12.5*(Vm+0.044));     % potassium inactivation rate
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    m_inf = alpha_m./(alpha_m+beta_m);      % sodium activation variable
    
    tau_h = 1./(alpha_h+beta_h);            % sodium inactivation time constant
    h_inf = alpha_h./(alpha_h+beta_h);      % sodium inactivation variable
    
    tau_n = 1./(alpha_n+beta_n);            % potassium activation time constant
    n_inf = alpha_n./(alpha_n+beta_n);      % potassium activation variable
    
    mca_inf = 1./(1+exp(-(Vm+0.052)/0.0074));   % Ca T-current activation variable
    
    hca_inf = 1./(1+exp(500*(Vm+0.076)));       % Ca T-current inactivation variable
    
    % Below the Ca T-current time constant has two terms, one for Vm<-80mV,
    % the other for Vm > -80mV. The expressions in parenthesis, such as
    % (Vm<-0.080) yield 1 or 0 depending on the value of Vm, so either the
    % first line is used or the second term line used.
    tau_hca = 1e-3*exp(15*(Vm+0.467)).*(Vm<-0.080) ...
        + 1e-3*(28+exp(-(Vm+0.022)/0.0105)).*(Vm>=-0.080);
    
    m(i) = m_inf;    % Update m, assuming time constant is neglible.
    
    h(i) = h(i-1) + dt*(h_inf - h(i-1))/tau_h;    % Update h
    
    n(i) = n(i-1) + dt*(n_inf - n(i-1))/tau_n;    % Update n
    
    mca(i) = mca_inf;    % Update mca, assuming time constant is negligible
    hca(i) = hca(i-1) + dt*(hca_inf - hca(i-1))/tau_hca;    % Update hca
    
    % Now calculate instantaneous values of the conductances
    G_Na_now = G_Na*m(i)*m(i)*m(i)*h(i);
    
    G_K_now = G_K*n(i)*n(i)*n(i)*n(i);
    
    G_CaT_now = G_CaT*mca(i)*mca(i)*hca(i);
    
    % Now update the membrane potential
    V(i) = Vm + dt*(G_L*(E_L-Vm) + G_Na_now*(E_Na-Vm) ...
        +G_K_now*(E_K-Vm) + G_CaT_now*(E_Ca-Vm) + I(i-1) )/Cm;
end
