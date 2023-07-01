% Figure_5_11A.m
% This code runs through multiple trials of applied current, commencing
% each trial from the final state of the prior trial.
%
% This model is the Hodgkin-Huxley model in new units.
%
% This code is used to produce Figure 5.11A of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dt = 2e-8;          % time-step for integration (sec)


%% Neuron parameters

V_L = -0.060;       % leak reversal potential (V)
E_Na = 0.045;       % reversal for sodium channels (V)
E_K = -0.082;       % reversal for potassium channels (V)
V0 = -0.065;

G_L = 30e-9;        % specific leak conductance (S)
G_Na = 12e-6;       % specific sodium conductance (S)
G_K = 3.6e-6;         % specific potassium conductance (S)

Cm = 100e-12;       % specific membrane capacitance (F)

Ibase = 0.7e-9;
tmax=0.25;             % maximum time of simulation (s)
t=0:dt:tmax;        % time vector
V=zeros(size(t));           % membrane potential vector
n=zeros(size(t));       % n: potassium activation gating variable
m=zeros(size(t));       % m: sodium activation gating variable
h=zeros(size(t));       % h: sodim inactivation gating variable
V(1) = V_L;
h(1) = 0.5;

Itot=zeros(size(t));    % in case we want to plot and look at the total current
I_Na=zeros(size(t));    % record sodium curret
I_K=zeros(size(t));     % record potassium current
I_L=zeros(size(t));     % record leak current

%% Now add current pulses at different points on the cycle and analyze the
%  change in response due to the pulse.

Ipulse_amp = 10e-12;

Npulses = 200;
Ipulse = 0;
shift = zeros(1,Npulses);

for trial = 0:Npulses;  % INitially trial 0 is used to get the baseline
    
    if ( trial == 1 )
        i_tpulse = i_startphase + floor((i_stopphase-i_startphase)*[1:Npulses]/Npulses);
    end
    
    Iapp=Ibase*ones(size(t)); % Applied current, relevant in current-clamp mode
    if ( trial > 0 )
        Iapp(i_tpulse(trial):i_tpulse(trial)+round(0.005/dt)) = ...
            Iapp(i_tpulse(trial):i_tpulse(trial)+round(0.005/dt)) + Ipulse_amp;
    end
    lastspiketime = 0;
    inspike = 0;
    
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

    %% Detect the initiation time of individual bursts
    inspike = 0;
    tspike = [];
    Nspikes = 0;
    for i = 1:length(t)
        if ( inspike == 0 ) && ( V(i) > -0.010 )
            inspike = 1;
            tspike(end+1) = t(i);
            Nspikes = Nspikes + 1;
        end
        if (inspike == 1 ) && ( V(i) < -0.05 )
            inspike = 0;
        end
    end
    %% Now decide where a phase of "zero" corresponds to in the oscillation
    %  and calculate the time where this occurs -- time should be well after any
    %  initial transients.
    if ( trial == 0 )
        i_startphase = round(tspike(3)/dt);
        i_stopphase = round(tspike(4)/dt);
        period = (tspike(4) - tspike(3));
    else
        shift(trial) = 2*pi*(1 - (tspike(4) - tspike(3))/period );
        if ( shift(trial) < -pi )
            shift(trial) = shift(trial) + 2*pi;
        end
        if ( shift(trial) > pi )
            shift(trial) = shift(trial) - 2*pi;
        end
        
        data = [shift(trial) Nspikes tspike(4)]
    end
    
    
    
    %% Now plot the graphs
    if ( mod(trial,10) == 0 )
        set(0,'DefaultLineLineWidth',2,...
            'DefaultLineMarkerSize',8, ...
            'DefaultAxesLineWidth',2, ...
            'DefaultAxesFontSize',14,...
            'DefaultAxesFontWeight','Bold');
        
        figure(trial+1)
        clf
        plot(t,V,'k');
        xlabel('Time (sec)')
        ylabel('Membrane potential (V)')
        axis([tspike(3)-0.005 tspike(4)+0.005 -0.085 0.050])
        drawnow
    end
end

figure(5)
clf
plot([1:Npulses]*2*pi/Npulses,shift,'k')
hold on
plot([1:Npulses]*2*pi/Npulses,[1:Npulses]*0,'k:')
xlabel('Phase of pulse')
ylabel('Phase shift')
axis([0 2*pi -0.15 0.18])
set(gca,'YTick',[-0.1 0 0.1 0.2])
set(gca,'XTick',[0 pi 2*pi])
set(gca,'XTickLabel',{'0' '\pi' '2\pi'})
