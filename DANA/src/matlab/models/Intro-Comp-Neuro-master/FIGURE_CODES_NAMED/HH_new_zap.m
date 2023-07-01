% HH_new_zap.m
% This code applied an oscillating current of increasing frequency to the
% Hodgkin-Huxley model neuron.
% Type-II behavior of the Hodgkin-Huxley model is revealed by a resonance
% in the resulting "ZAP" response.
%
% This model is the Hodgkin-Huxley model in modern SI units.
%
%
% This code is used to produce Figure 4.10 in the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dt = 2e-8;          % time-step for integration (sec)
tmax=2;             % maximum time of simulation (sec)
t=0:dt:tmax;        % time vector
Nt = length(t);     % number of time points

Imaxvec = [0.1e-9 0.2e-9];  % Amplitude of sine wave current
Ntrials = length(Imaxvec);  % Number of trials with different amplitude

fmin = 0;                   % Minimum (initial) oscillating frequency
fmax = 80;                  % Half of the maximum (final) oscillating frequency
df = (fmax-fmin)/(Nt-1);    % Half of the change in frequency per time-step
foft = fmin:df:fmax;        % Half of the frequency on a given time-step

%% Neuron parameters
V_L = -0.060;       % leak reversal potential (V)
E_Na = 0.045;       % reversal for sodium channels (V)
E_K = -0.082;       % reversal for potassium channels (V)
V0 = -0.065;        % initial membrane potential

G_L = 30e-9;        % leak conductance (S)
G_Na = 12e-6;       % sodium conductance (S)
G_K = 3.6e-6;       % potassium conductance (S)

Cm = 100e-12;       % membrane capacitance (F)

V=zeros(size(t));       % membrane potential vector
n=zeros(size(t));       % n: potassium activation gating variable
m=zeros(size(t));       % m: sodium activation gating variable
h=zeros(size(t));       % h: sodim inactivation gating variable

%% Now loop through trials, each with a different amplitude
for trial = 1:Ntrials
    
    Imax = Imaxvec(trial)           % Amplitude of oscillation    
    Iapp = Imax*sin(2*pi*foft.*t);  % Oscillating applied current
    % Note that the simple method of producing an oscillating current with
    % linearly varying frequency means the time difference between peaks is
    % 1/(2<f>) where <f> is the mean value of foft between the two peaks
    
    %% General initialization
    V(1) = -0.065;          % set the inititial value of voltage
    m(1) = 0.05;            % initialize sodium activation
    h(1) = 0.6;             % initalize sodium inactivation
    n(1) = 0.31;            % initialize potassium activation
    Itot=zeros(size(t));    % record the total current
    I_Na=zeros(size(t));    % record sodium curret
    I_K=zeros(size(t));     % record potassium current
    I_L=zeros(size(t));     % record leak current
    
    for i = 2:length(t);    % now see how things change through time
        
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
        
        % Update gating variables using the Forward Euler method
        m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
        h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h        
        n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
        
        % Calculate currents from gating variables
        I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i-1)); % total sodium current        
        I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i-1)); % total potassium current        
        I_L(i) = G_L*(V_L-V(i-1));    % Leak current is straightforward        
        Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i); % total current is sum of leak + active channels + applied current

        % Update membrane potential
        V(i) = V(i-1) + Itot(i)*dt/Cm;        % Update the membrane potential, V.
        
    end
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    % figure(1) becomes Figure 4.10 in the book
    figure(1)
    
    % First plot the applied current profile, but just once since it is scaled
    if ( trial == 1)
        clf
        subplot('Position',[0.12 0.72 0.84 0.23])
        plot(t(100:100:end)/1000,Iapp(100:100:end)/Imax,'k');
        ylabel('I_{app}')
        set(gca,'YTick',[-1 0 1])
        set(gca,'YTickLabel',{'-I_{max}', '0', 'I_{max}'})
    end
    
    % The second and third panels contain the membrane potential versus
    % time
    subplot('Position',[0.12 0.72-0.29*trial 0.84 0.23])
    plot(t(100:100:end),1000*V(100:100:end),'k');
    ylabel('Vm (mV)')
    legstring = strcat('I_{max} = ',num2str(Imax*1e9),'nA')
    legend(legstring)
    if ( trial == 1 )
        axis([0 tmax -74 -66]) % zoom in for subthreshold response
    else
        axis([0 tmax -85 35])
    end
end

xlabel('Time (sec)')
