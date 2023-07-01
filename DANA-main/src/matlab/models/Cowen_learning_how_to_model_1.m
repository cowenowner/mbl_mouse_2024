%% Get a sense of what a negative exponential is - as it is important for
% understanding models of decay - so common in neuro and chemistry.
% CHAPTER 2
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 0.1:.001:10;

figure
plot(x,exp(-(x))*1)
hold on;plot(x,exp(-(x))*.4)
hold on;plot(x, exp(-(x))*1.4)
% the equation for voltage change - all passive properties.
dt = 0.001;
t_s = 0:dt:0.05;
E_L = -0.070;       % Leak potential (steady state or resting potential)
V_V = zeros(size(t_s));
V_V(1) = -.050;
tau_s = 0.010;
figure
plot(t_s,E_L + (V_V(1)-E_L)*exp(-t_s/tau_s))
xlabel('s');ylabel('V')
%% Do the same but with differential.
% C_m = .2;
V_V = zeros(size(t_s));
V_V(1) = -.050;

for ii = 2:length(t_s)
    V_V(ii) = V_V(ii-1) + dt*(E_L-V_V(ii-1))/tau_s;
end
figure
plot(t_s,V_V)
xlabel('s');ylabel('V')

%% Exercise 2.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leaky integrate and fire model...
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see LIF_model.m
% 2.1a
dt = 0.001;
t_max = 2;
tau = 0.010;                % membrane time constant
t_s = 0:dt:t_max;
Cm_nF = 2e-9;               % total membrane capacitance nano Farads
G_L = Cm/tau;               % total membrane conductance (leak conductance)
E_L = -0.07;
R_m_MO = 5;
V_V = zeros(size(t_s));
V_reset = -0.065;
V_V(1) = E_L;
I_app = zeros(size(t_s)) + 12e-12;
I_app = zeros(size(t_s)) + 512e-12;
I_app = zeros(size(t_s)) + 2.1e-10;

th_V = -0.05;
r = randn(size(t_s))*1810e-12 ;
s = zeros(size(t_s));

for ii = 2:length(t_s)
%     V_V(ii) = V_V(ii-1) + dt*(E_L-V_V(ii-1))/tau_s + r(ii);
    V_V(ii) = V_V(ii-1) + dt*( I_app(ii) + G_L*(E_L-V_V(ii-1)))/Cm_nF;
%             V(i) = V(i-1) + dt*(I(i) +G_L*(E_L-V(i-1)))/Cm;

    if V_V(ii) > th_V
        V_V(ii) = V_reset;
        s(ii) = 1;
    end
    if V_V(ii) < V_reset*2
        V_V(ii) = V_reset;
    end
end
figure
plot(t_s,V_V)
plot_markers(t_s(s>0))

xlabel('s');ylabel('V')
%% 1b Determine min applied current to get AP...
% Ith = GL(Vth-EL)
Ith = G_L*(th_V - E_L)
% 1c = just goes through Iapp - won't bother.
% 1d Computer the FI curve.
Iapp = linspace(Ith*.2,Ith*2,100);
ISI_s = tau*log(Iapp*R_m_MO+E_L-V_reset) - tau*log(Iapp*R_m_MO+E_L-th_V);
figure
plot(Iapp,1./ISI_s); axis tight
xlabel('iApp');ylabel('Hz')
% The above is wrong somehow. Negative values for rate and ISI - cant' be
% rihgt.

% Tutorial 2.2 Refractory period.


