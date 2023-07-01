% coupledPIR.m 
% A model of a half-center oscillator by coupling together two neurons with
% inhibitory synaptic connections. 
%
% Each model neuron contains a T-type Calcium current to generate a
% post-inhibitory rebound as a models of thalamic relay cells.
%
% The model is based on the paper by Wang and Rinzel.
%
% This code was used to produce Figure 5.9A in the textbook
% An Introductory Course in Computational Neuroscience 
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This two-cell model is written without using two-row matrices for variables, 
% which means that each cell has a separate variable name for each of its
% variables. This means the code is lengthier than necessary, but does not
% impact its memory or runtime.
%
% Parameters have a single name, as they are identical across cells.
%
clear;              % clear all variables and parameters from memory
dt = 0.0001;        % time step in sec
tmax=1;             % total simulation time in sec
t=0:dt:tmax;        % vector of time points

istart = 0.25;          % time applied current starts (sec)
ilength=tmax-istart;    % length of applied current pulse (sec)
I01= -5e-9;         % Base current to cell 1 before pulse
Ie1= 0e-9;          % magnitude of applied current pulse
I02= -8e-9;         % Base current to cell 2 before pulse
Ie2= 0e-9;          % magnitude of applied current pulse

E_L = -0.070;   % leak reversal potential
E_Na = 0.055;   % reversal for sodium channels
E_K = -0.090;   % reversal for potassium channels
E_Ca = 0.120    % reversal potential for Ca current
E_Syn = -0.080  % reversal potential of inhibitory synapse
Vspike = 0.0    % voltage at which synapse is activated

g_L = 1.0e-6;     % specific leak conductance in Siemens per mm-square
g_Na = 0.32e-3;   % specific sodium conductance
g_K = 0.18e-3;    % specific potassium conductance
g_CaT = 20e-6;    % T-type calcium conductance
g_syn12 = 20e-6   % synaptic conductance from cell 1 to 2
g_syn21 = 20e-6   % synaptic conductance from cell 2 to 1
tau_syn = 0.005   % synaptic time constant
cm = 10e-9;     % ecific membrane capacitance in Farads per mm-square

%% Begin with all variables for cell 1
% Define the arrays
V1=zeros(size(t));          % membrane potential
V1(1) = E_L;                % set the initial value of voltage     
I1_L= zeros(size(t));       % leak current
I1_Na= zeros(size(t));      % sodium current
I1_K= zeros(size(t));       % potassium current
I1_CaT = zeros(size(t));    % T-type calcium current

% Gating variables for cell 1
n1=zeros(size(t));      % n: potassium activation gating variable
n1(1) = 0.0;            % start off at zero
m1=zeros(size(t));      % m: sodium activation gating variable
m1(1) = 0.0;            % start off at zero
h1=zeros(size(t));      % h: sodium inactivation gating variplot(t,V)able
h1(1) = 0.5;            % start off at zero
mca1=zeros(size(t));    % CaT current activation gating variable
mca1(1) = 0.0;          % start at zero
hca1=zeros(size(t));    % CaT current inactivation gating variable
hca1(1) = 0.0;          % start at zero

% Set up applied current to cell 1
I1app=zeros(size(t));   
for i = 1:round(istart/dt)
    I1app(i) = I01;     % baseline (negative here)
end
for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
    I1app(i) = Ie1;     % applied current (in fact zero)
end
for i = round((istart+ilength)/dt):length(I1app)
    I1app(i) = I01;     % this is not used
end

I1tot=zeros(size(t)); % in case we want to plot and look at the total current

%% Repeat with variables for cell 2
% Define the arrays
V2=zeros(size(t));          % membrane potential for cell 2
V2(1) = E_L;                % set the inititial value of voltage     
I2_L= zeros(size(t));       % leak current
I2_Na= zeros(size(t));      % sodium current
I2_K= zeros(size(t));       % potassium current
I2_CaT = zeros(size(t));    % T-type calcium current

% Gating variables for cell 2
n2=zeros(size(t));      % n: potassium activation gating variable
n2(1) = 0.0;            % start off at zero
m2=zeros(size(t));      % m: sodium activation gating variable
m2(1) = 0.0;            % start off at zero
h2=zeros(size(t));      % h: sodim inactivation gating variplot(t,V)able
h2(1) = 0.0;            % start off at zero

mca2=zeros(size(t));    % CaT current activation gating variable
mca2(1) = 0.0;          % initialize at zero
hca2=zeros(size(t));    % CaT current inactivation gating variable
hca2(1) = 0.0;          % initialize at zero

% Set up applied current to cell 2
I2app=zeros(size(t));   % Applied current vector
for i = 1:round(istart/dt)
    I2app(i) = I02;
end
for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
    I2app(i) = Ie2;
end
for i = round((istart+ilength)/dt):length(I2app)
    I2app(i) = I02;
end
    
I2tot=zeros(size(t));   % to store the total current

syn1=zeros(size(t));    % synaptic gating variable from spikes in cell 1
syn2=zeros(size(t));    % synaptic gating variable from spikes in cell 2
spikes1 = zeros(size(t));   % store spikes of cell 1
spikes2 = zeros(size(t));   % store spikes of cell 2
spike1now = 0;              % set to 1 when cell 1 is in a spike
spike2now = 0;              % set to 1 when cell 2 is in a spike

% nspikewidth is number of time-points to update conductance following a spike
nspikewidth = round(tau_syn*10/dt);  

%% Now loop through time and simulate the two cells
for i = 2:length(t); % now see how things change through time

    %% Calculate currents and voltage for cell 1 first
    I1_L(i-1) = g_L*(E_L-V1(i-1));  % Leak current
    
    % Notice that synaptic current to cell 1 depends on syn2, which is
    % impacted by spikes from cell 2, but also on V1, membrane potential of
    % cell 1.
    I1syn(i-1) = g_syn21*syn2(i-1)*(E_Syn-V1(i-1));
    
    Vm = V1(i-1)*1000; % converts voltages to mV as needed in the equations on p.224 of Dayan/Abbott
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta. Written out from Dayan/Abbott, units are 1/ms.
    if ( Vm == -35 ) 
        alpha_m = 1;
    else 
        alpha_m = 0.1*(Vm+35)/(1-exp(-0.1*(Vm+35)));
    end
    beta_m = 4*exp(-(Vm+60)/18);
    alpha_h = 0.35*exp(-0.05*(Vm+58));
    beta_h = 5/(1+exp(-0.1*(Vm+28)));
    
    if ( Vm == -34 ) 
       alpha_n = 0.05/0.1;
    else
        alpha_n = 0.05*(Vm+34)/(1-exp(-0.1*(Vm+34)));
    end
    beta_n = 0.625*exp(-0.0125*(Vm+44));
     
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    m_inf = alpha_m/(alpha_m+beta_m);   % immediate change in m (no time constant)
    
    tau_h = 1e-3/(alpha_h+beta_h);      % time constant converted from ms to sec
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1e-3/(alpha_n+beta_n);      % time constant converted from ms to sec
    n_inf = alpha_n/(alpha_n+beta_n);   
    
    % for the Ca_T current gating variables are given by formulae for the 
    % steady states and time constants:
    mca_inf = 1/(1+exp(-(Vm+52)/7.4));
    hca_inf = 1/(1+exp((Vm+76)/2));
    if ( Vm < -80 ) 
        tau_hca = 1e-3*exp((Vm+467)/66.6);
    else
        tau_hca = 24e-3+1e-3*119/(1+exp((Vm+70)/3));
    end
    
    % Now update all gating variables using "steady states" and time
    % constants
    m1(i) = m_inf;    % Update m, assuming time constant is neglible.
    h1(i) = h_inf - (h_inf-h1(i-1))*exp(-dt/tau_h);    % Update h
    n1(i) = n_inf - (n_inf-n1(i-1))*exp(-dt/tau_n);    % Update n
    mca1(i) = mca_inf;
    hca1(i) = hca_inf - (hca_inf-hca1(i-1))*exp(-dt/tau_hca);
    
    % Update conductances then currents using gating variables
    g_Na_now = g_Na*m1(i)*m1(i)*m1(i)*h1(i);
    I1_Na(i-1) = g_Na_now*(E_Na-V1(i-1)); % total sodium current
    
    g_K_now = g_K*n1(i)*n1(i)*n1(i)*n1(i);
    I1_K(i-1) = g_K_now*(E_K-V1(i-1)); % total potassium current
    
    g_CaT_now = g_CaT*mca1(i)*mca1(i)*hca1(i);
    I_CaT(i-1) = g_CaT_now*(E_Ca-V1(i-1)); % Calcium T-type current
    
    % Total current to cell 1
    I1tot(i-1) = I1_L(i-1)+I1_Na(i-1)+I1_K(i-1) ... % total current is sum of leak + active  
                +I1_CaT(i-1) + I1syn(i-1) +I1app(i-1);           % + synaptic + applied current
     
    g_Tot = g_L+g_Na_now+g_K_now+g_CaT_now+g_syn21*syn2(i-1);
    V_inf = (g_L*E_L + g_Na_now*E_Na + g_K_now*E_K  + ...
            g_CaT_now*E_Ca + g_syn21*syn2(i-1)*E_Syn + I1app(i-1))/g_Tot;
    
    % Update membrane potential of cell 1 using the stable exponential
    % method
    V1(i) = V_inf - (V_inf-V1(i-1))*exp(-dt*g_Tot/cm);  % Update the membrane potential, V.

    %% Calculate currents and voltage for cell 2 next
    I2_L(i-1) = g_L*(E_L-V2(i-1));
    I2syn(i-1) = g_syn12*syn1(i-1)*(E_Syn-V2(i-1));
    
    Vm = V2(i-1)*1000; % converts voltages to mV as needed in the equations on p.224 of Dayan/Abbott
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta. Written out from Dayan/Abbott, units are 1/ms.  
    
    if ( Vm == -35 ) 
        alpha_m = 1;
    else 
        alpha_m = 0.1*(Vm+35)/(1-exp(-0.1*(Vm+35)));
    end
    beta_m = 4*exp(-(Vm+60)/18);

    alpha_h = 0.35*exp(-0.05*(Vm+58));
    beta_h = 5/(1+exp(-0.1*(Vm+28)));
    
    if ( Vm == -34 ) 
       alpha_n = 0.05/0.1;
    else
        alpha_n = 0.05*(Vm+34)/(1-exp(-0.1*(Vm+34)));
    end
    beta_n = 0.625*exp(-0.0125*(Vm+44));
     
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1e-3/(alpha_h+beta_h);      % time constant converted from ms to sec
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1e-3/(alpha_n+beta_n);      % time constant converted from ms to sec
    n_inf = alpha_n/(alpha_n+beta_n);   
    
    % for the Ca_T current gating variables are given by formulae for the 
    % steady states and time constants:
    
    mca_inf = 1/(1+exp(-(Vm+52)/7.4));

    hca_inf = 1/(1+exp((Vm+76)/2));
    if ( Vm < -80 ) 
        tau_hca = 1e-3*exp((Vm+467)/66.6);
    else
        tau_hca = 24e-3+1e-3*119/(1+exp((Vm+70)/3));
    end
    
    m2(i) = m_inf;    % Update m, assuming time constant is neglible.
    
    h2(i) = h_inf - (h_inf-h2(i-1))*exp(-dt/tau_h);    % Update h
    
    n2(i) = n_inf - (n_inf-n2(i-1))*exp(-dt/tau_n);    % Update n
    
    mca2(i) = mca_inf;
    hca2(i) = hca_inf - (hca_inf-hca2(i-1))*exp(-dt/tau_hca);
    
    g_Na_now = g_Na*m2(i)*m2(i)*m2(i)*h2(i);
    I2_Na(i-1) = g_Na_now*(E_Na-V2(i-1)); % total sodium current
    
    g_K_now = g_K*n2(i)*n2(i)*n2(i)*n2(i);
    I2_K(i-1) = g_K_now*(E_K-V2(i-1)); % total potassium current
    
    g_CaT_now = g_CaT*mca2(i)*mca2(i)*hca2(i);
    I2_CaT(i-1) = g_CaT_now*(E_Ca-V2(i-1)); % Calcium T-type current
    
    
    I2tot(i-1) = I2_L(i-1)+I2_Na(i-1)+I2_K(i-1) ...% total current is sum of synaptic 
        +I2_CaT(i-1) + I2syn(i-1) + I2app(i-1);  % + leak + active channels + applied current
     
    g_Tot = g_L+g_Na_now+g_K_now+g_CaT_now+g_syn12*syn1(i-1);
    V_inf = (g_L*E_L + g_Na_now*E_Na + g_K_now*E_K  + ...
            g_CaT_now*E_Ca + g_syn12*syn1(i-1)*E_Syn + I2app(i-1))/g_Tot;
                   
    V2(i) = V_inf - (V_inf-V2(i-1))*exp(-dt*g_Tot/cm);  % Update the membrane potential, V.

    %% Now check for spikes and update synaptic gating variables
    % First for cell 1
    if ( V1(i) > Vspike ) && (spike1now == 0 )  % if V1 is at a spiking level, but the spike is not already detected
        spike1now = 1;                          % detect this spike
        spikes1(i) = 1;                         % record this spike time
        for j = 0:nspikewidth                   % for many synaptic time constants in the future
            if i+j <= length(syn1)              % so long as we are still inside the total simulation time
                syn1(i+j) = syn1(i+j) + dt*j/tau_syn*exp(-dt*j/tau_syn);    % extra conductance from this spike
            end
        end
    end
    if ( V1(i) < Vspike - 0.010 )               % once the spike is over
        spike1now = 0;                          % set this to zero so we are ready for the next spike.
    end
    
    % Now do the same for spikes from neuron 2.
    if ( V2(i) > Vspike ) && (spike2now == 0 ) 
        spike2now = 1;
        spikes2(i) = 1;
        for j = 0:nspikewidth
            if i+j <= length(syn1) 
                syn2(i+j) = syn2(i+j) + dt*j/tau_syn*exp(-dt*j/tau_syn);
            end
        end
    end
    if ( V2(i) < Vspike - 0.010 ) 
        spike2now = 0;
    end
    
end

%% Finally plot membrane potentials as a function of time
figure(2)
% Membrane potential of cell 1 in the upper panel
subplot(2,1,1)      % upper panel
plot(t,V1,'k')
xlabel('Time (sec)')
ylabel('Cell 1 Vm')
hold on
% Membrane potential of cell 2 in the lower panel
subplot(2,1,2)      % lower panel
plot(t,V2,'k')
xlabel('Time (sec)')
ylabel('Cell 2 Vm')




    