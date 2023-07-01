% PIR_Vdependence.m
% This model contains a T-type Calcium current to generate a
% post-inhibitory rebound as a models of thalamic relay cells.
%
%   The code provides an implementation of the gating variables for 
%   Tutorial 4.2 of Chapter 4 in the textbook,
%   An Introductory Course in Computational Neuroscience 
%   by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

E_L = -0.070;       % leak reversal potential
E_Na = 0.055;       % reversal for sodium channels
E_K = -0.090;       % reversal for potassium channels
E_Ca = 0.120        % reversal potential for Ca current

G_L = 10e-9;        % leak conductance in Siemens 
G_Na = 3.6e-6;      % sodium conductance in Siemens
G_K = 1.6e-6;       % potassium conductance in Siemens
G_CaT = 0.22e-6;    % T-type calcium conductance in Siemens

Cm = 0.1e-9;        % membrane capacitance in Farads 
Vm = -0.1:0.001:0;  % range of membrane potential (the variable here)

% alpha_m is sodium activation rate constant
if ( Vm == -35 )                            % to prevent division by zero
    alpha_m = 1e3;
else
    alpha_m = 1e5*(Vm+0.035)./(1-exp(-100*(Vm+0.035)));
end
beta_m = 4000*exp(-(Vm+0.060)./0.018);      % sodium deactivation rate

alpha_h = 350*exp(-50*(Vm+0.058));          % sodium inactivation rate
beta_h = 5000./(1+exp(-100*(Vm+0.028)));    % sodium deinactivation rate

if ( Vm == -34 )                            % prevent division by zero
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

%% Now plot the functions after setting up graphics parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf
plot(Vm,mca_inf)
hold on
plot(Vm,hca_inf)

figure(2)
plot(Vm,tau_hca)



