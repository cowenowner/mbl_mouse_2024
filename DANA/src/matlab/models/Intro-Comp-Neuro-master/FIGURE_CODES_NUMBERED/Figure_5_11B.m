% Figure_5_11B.m
%
% This model is the Connors-Stevens model, similar to Hodgkin-Huxley, but
% more like neurons in the cortex, being type-I.
% See Dayan and Abbott Sect 5.5, pp 166-172 then Sect 6.1-2, pp.196-198 and Sect 6.6 p.224.
%
% This code is used to produce Figure 5.11B of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dt = 2e-6;
tmin = 0;
tmax= 2;

t=tmin:dt:tmax;    % time vector

istart = 0; % time applied current starts
ilength=tmax;   % length of applied current pulse
nstart = floor((istart-tmin)/dt)+1;
nstop = floor((istart+ilength-tmin)/dt)+1;
Ibase = 0.82e-9;    % applied current

V_L = -0.017;   % leak reversal potential
E_Na = 0.055;   % reversal for sodium channels
E_K = -0.072;   % reversal for potassium channels
E_A = -0.075;   % reversal for A-type current

G_L = 3e-8;     % leak conductance
G_Na = 1.2e-5;  % sodium conductance
G_K = 2e-6;     % potassium conductance
G_A = 4.77e-6;  % A-tpe potassium conductance

Cm = 0.1e-9;     % membrane capacitance

V=zeros(size(t)); % voltage vector
V(1) = V_L;    % set the inititial value of voltage

n=zeros(size(t));   % n: potassium activation gating variable
n(1) = 0.0;         % start off at zero
m=zeros(size(t));   % m: sodium activation gating variable
m(1) = 0.0;         % start off at zero
h=zeros(size(t));   % h: sodim inactivation gating variable
h(1) = 0.0;         % start off at zero

a=zeros(size(t));   % A-current activation gating variable
a(1) = 0.0;         % start off at zero
b=zeros(size(t));   % A-current inactivation gating variable
b(1) = 0.0;         % start off at zero

Itot=zeros(size(t)); % in case we want to plot and look at the total current
I_Na=zeros(size(t));
I_K=zeros(size(t));
I_A=zeros(size(t));
I_L=zeros(size(t));
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
    
    for i = 2:length(t); % now see how things change through time
        
        Vm = V(i-1); % converts voltages to mV as needed in the equations on p.224 of Dayan/Abbott
        
        % Sodium and potassium gating variables are defined by the
        % voltage-dependent transition rates between states, labeled alpha and
        % beta. Written out from Dayan/Abbott, units are 1/ms.
        
        alpha_m = 3.80e5*(Vm+0.0297)/(1-exp(-100*(Vm+0.0297)));
        beta_m = 1.52e4*exp(-55.6*(Vm+0.0547));
        
        alpha_h = 266*exp(-50*(Vm+0.048));
        beta_h = 3800/(1+exp(-100*(Vm+0.018)));
        
        alpha_n = 2e4*(Vm+0.0457)/(1-exp(-100*(Vm+0.0457)));
        beta_n = 250*exp(-12.5*(Vm+0.0557));
        
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        
        tau_m = 1/(alpha_m+beta_m);      % time constant converted from ms to sec
        m_inf = alpha_m/(alpha_m+beta_m);
        
        tau_h = 1/(alpha_h+beta_h);      % time constant converted from ms to sec
        h_inf = alpha_h/(alpha_h+beta_h);
        
        tau_n = 1/(alpha_n+beta_n);      % time constant converted from ms to sec
        n_inf = alpha_n/(alpha_n+beta_n);
        
        m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
        h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
        n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        
        % For the A-type current gating variables, instead of using alpha and
        % beta, we just use the steady-state values a_inf and b_inf along with
        % the time constants tau_a and tau_b that are found empirically
        a_inf = (0.0761*exp(31.4*(Vm+0.09422))/(1+exp(34.6*(Vm+0.00117))))^(1/3.0);
        tau_a = 0.3632e-3 + 1.158e-3/(1+exp(49.7*(Vm+0.05596)));
        
        b_inf = (1/(1+exp(68.8*(Vm+0.0533))))^4;
        tau_b = 1.24e-3 + 2.678e-3/(1+exp(62.4*(Vm+0.050)));
        
        a(i) = a(i-1) + (a_inf-a(i-1))*dt/tau_a;    % Update a
        b(i) = b(i-1) + (b_inf-b(i-1))*dt/tau_b;    % Update b
        
        I_L(i) = G_L*(V_L-V(i-1));
        
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i-1)); % total sodium current
        
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i-1)); % total potassium current
        
        I_A(i) = G_A*a(i)*a(i)*a(i)*b(i)*(E_A-V(i-1)); % total A-type current
        
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+I_A(i)+Iapp(i); % total current is sum of leak + active channels + applied current
        
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
        axis([0.65 1.05 -0.075 0.05])
        drawnow
    end
end

figure()
plot([1:Npulses]*2*pi/Npulses,shift,'k')
hold on
plot([1:Npulses]*2*pi/Npulses,[1:Npulses]*0,'k:')
xlabel('Phase of pulse')
ylabel('Phase shift')
axis([0 2*pi -0.1 0.2])
set(gca,'YTick',[-0.1 0 0.1 0.2])
set(gca,'XTick',[0 pi 2*pi])
set(gca,'XTickLabel',{'0' '\pi' '2\pi'})


