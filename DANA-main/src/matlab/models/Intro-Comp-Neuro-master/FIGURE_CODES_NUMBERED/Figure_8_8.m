% Figure_8_8.m
% This model is the Connors-Stevens model, similar to Hodgkin-Huxley, but
% more like neurons in the cortex, being type-I.
% See Dayan and Abbott Sect 5.5, pp 166-172 then Sect 6.1-2, pp.196-198 and Sect 6.6 p.224.
%
% In this version, synaptic inputs arrive separately then synchronously. 
% With the default parameters the model neuron is a coincidence detector,
% only responding to the combined synchronous input.
% In separate trials, the synaptic strength is altered from its default
% value, or the A-type potassium conductance is altered from its default
% value. These modifications impair the neuron's behavior.
%
%  This code is used to produce Figure 8.8 in Chapter 8 of the textbook
%  An Introductory Course in Computational Neuroscience
%  by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dt = 0.000005;      % time-step for simulations
tmax=1;             % maximum time to simulate each trial
t=0:dt:tmax;        % time vector
Nt = length(t);     % Number of time-steps

V_L = -0.017;   % leak reversal potential
E_Na = 0.055;   % reversal for sodium channels
E_K = -0.072;   % reversal for potassium channels
E_A = -0.075;   % reversal for A-type current
E_syn = 0.000;  % synaptic reversal potential (excitatory)

G_L = 3e-8;     % specific leak conductance
G_Na = 1.2e-5;  % specific sodium conductance
G_K = 2e-6;     % specific potassium conductance
G_A0 = 5e-6;    % specific A-tpe potassium conductance (default value)
G_syn_max0 = 1e-7;  % synaptic input conductance

Cm = 0.1e-9;    % specific membrane capacitance

taus = 0.002;   % synaptic time constant of inputs

t=0:dt:tmax;    % time vector
Nt = length(t); % Number of time-steps

%% Produce trial-specific values of parameters
Ntrials = 6;                            % Number of trials
% Initialize the trial-specific values with the default values
G_A_trials = G_A0*ones(1,Ntrials);     
G_syn_max_trials = G_syn_max0*ones(1,Ntrials);

G_syn_max_trials(1) = 5e-8;     % Reduce synaptic conductance on trial 1
G_syn_max_trials(3) = 1.5e-7;   % Increase synaptic conductance on trial 3

G_A_trials(4) = 7.5e-6;         % Increase A-type conductance on trial 4
G_A_trials(6) = 2.5e-6;         % Reduce A-type conductance on trial 6

%% Set up the presynaptic spike times to be used, then the synaptic input
inspikes1 = zeros(1,Nt);        % Initialize to zero input spikes
inspikes2 = zeros(1,Nt);        % Initialize to zero input spikes
inspikes1(floor(Nt/4)) = 1;     % Spike time for input 1
inspikes1(floor(3*Nt/4)) = 1;   % Second (synchronous) spike for input 1
inspikes2(floor(Nt/2)) = 1;     % Spike time for input 2
inspikes2(floor(3*Nt/4)) = 1;   % Second (synchronous) spike for input 2

s1 = zeros(1,Nt);               % Vector for synaptic gating variable of input 1
s2 = zeros(1,Nt);               % Vector for synaptic gating variable of input 2
for i = 2:Nt
    s1(i) = s1(i-1)*(1-dt/taus);    % Decays with time constant taus
    s1(i) = s1(i) + inspikes1(i);   % Increments at spike times of input 1
    s2(i) = s2(i-1)*(1-dt/taus);    % Decays with time constant taus
    s2(i) = s2(i) + inspikes2(i);   % Increments at spike times of input 2
end

%% Plot the inputs s1 above s2 at the top of each column of figures
figure(1)
clf
subplot('Position',[0.12 0.9 0.36 0.07])
plot(t,s1,'k')
set(gca,'YTick',[0 1])
ylabel('Input 1')

subplot('Position',[0.12 0.76 0.36 0.07])
plot(t,s2,'k')
set(gca,'YTick',[0 1])
ylabel('Input 2')

subplot('Position',[0.61 0.9 0.36 0.07])
plot(t,s1,'k')
set(gca,'YTick',[0 1])
ylabel('Input 1')

subplot('Position',[0.61 0.76 0.36 0.07])
plot(t,s2,'k')
set(gca,'YTick',[0 1])
ylabel('Input 2')

%% Now loop through trials with different values of parameters
for trial = 1:Ntrials
    
    G_syn_max = G_syn_max_trials(trial);    % Value of synaptic conductance for this trial
    G_A = G_A_trials(trial);                % Value of A-type conductance for this trial
        
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
    I_Na=zeros(size(t));    % Sodium current vector
    I_K=zeros(size(t));     % Delayed rectifier potassium current vector
    I_A=zeros(size(t));     % A-type potassium current vector
    I_L=zeros(size(t));     % Leak current vector
    
    for i = 2:length(t); % now see how things change through time
        
        Vm = V(i-1);    % Used to calculate gating variables from membrane potential
        
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
        
        I_syn(i) = G_syn_max*(s1(i)+s2(i))*(E_syn-V(i-1)); % total synaptic current
        
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+I_A(i)+I_syn(i); % total current is sum of leak + active channels + applied current
        
        V(i) = V(i-1) + Itot(i)*dt/Cm;        % Update the membrane potential, V.
                
    end
    
    %% Plot the membrane potential response for this trial    
    figure(1)
    col = ceil(trial/3);
    row = mod(trial-1,3)+1;
    subplot('Position',[col*0.49-0.37 0.76-row*0.22 0.36 0.15])
    plot(t,1000*V,'k')
    axis([0 1 -85 95])
    ylabel('V_{m} (mV)')
    if ( row == 3 ) 
        xlabel('Time (sec)')
    end
end

%% Label the figures A-G and add parameter values
annotation('textbox',[0.00 0.98 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A1')
annotation('textbox',[0.00 0.84 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A2')
annotation('textbox',[0.00 0.7 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.00 0.49 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.00 0.28 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','D')
annotation('textbox',[0.5 0.98 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A1')
annotation('textbox',[0.5 0.84 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A2')
annotation('textbox',[0.5 0.7 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','E')
annotation('textbox',[0.5 0.49 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','F')
annotation('textbox',[0.5 0.28 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','G')


annotation('textbox',[0.22 0.67 0.1 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','G_{syn}=50nS')
annotation('textbox',[0.22 0.45 0.1 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','G_{syn}=100nS')
annotation('textbox',[0.22 0.23 0.1 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','G_{syn}=150nS')
annotation('textbox',[0.72 0.67 0.1 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','G_{A}=7.5\muS')
annotation('textbox',[0.72 0.45 0.1 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','G_{A}=5\muS')
annotation('textbox',[0.72 0.23 0.1 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','G_{A}=2.5\muS')
