% graded_release.m
% Plots the functions for the steady-state conductance and the
% time-constant of a synapse with graded transmission as a function of
% presynaptic membrane potential.
% Then loads a presynaptic voltage trace from the Pinsky-Rinzel model in
% the file 'PR_VS.mat' and plots the synaptic conductance produced by such
% a bursting cell.
%
% This code is used to produce Figure 5.3 of the textbook,
% "An Introductory Course in Computational Neuroscience"
% by Paul Miller, Brandeis University, 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = -0.100:0.001:0.050;     % Range of presynaptic membrane potential

% The following are parameters from Prinz et al. Nat Neurosci 7:1345 (2004)
Vth = -0.035;
Vrange = 0.005;
tau_rise = 0.002;       % In Prinz et al tau_rise = 0
tau_decay = 0.025;
Gmax = 1e-9;
Gsyn_inf = Gmax./(1+exp(-(V-Vth)/Vrange));       % To plot m_h vs V
tau_Gsyn = tau_rise + (tau_decay-tau_rise)./(1 + exp((V-Vth)/Vrange)); % to plot tau_m_h vs V

figure(2)
clf
subplot('Position',[0.1 0.21 0.23 0.76])
plot(V*1000,Gsyn_inf*1e9,'k')
xlabel('Membrane potential (mV)')
ylabel('G_{\infty} (nS)')
axis([ -75 0 0 1])

subplot('Position',[0.43 0.21 0.23 0.76])
plot(V*1000,tau_Gsyn*1e3,'k')
xlabel('Membrane potential (mV)')
ylabel('\tau_{syn} (msec)')
axis([-75 0 0 26])

%% Now calculate the synaptic response to a bursting neuron

load('PR_VS.mat');  % Load somatic voltage trace from PR_euler_final.m

G = zeros(size(t)); % array to store postsynaptic conductance

for i = 2:length(t);    % step through time vector
    V = VS(i); % Presynaptic membrane potential (assume same as somatic VS)
    Gsyn_inf = Gmax./(1+exp(-(V-Vth)/Vrange));       % To plot m_h vs V
tau_Gsyn = tau_rise + (tau_decay-tau_rise)./(1 + exp((V-Vth)/Vrange)); % to plot tau_m_h vs V

    % Now update with the Forward Euler method
    G(i) = G(i-1) + dt*(Gsyn_inf-G(i-1))/tau_Gsyn;
end

subplot('Position',[0.76 0.21 0.23 0.76])
plot(t,G*1e9,'k')
hold on
plot(t,VS,'k:')
axis([0.850 0.925 -0.08 1])
xlabel('Time (sec)')
ylabel('G_{syn} (nS), V_{pre} (V)')
legend('G_{syn}', 'V_{pre}')

annotation('textbox',[0 1 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.35 1 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.67 1 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
