% CS_figure.m
%
% This model is the Connor-Stevens model, similar to Hodgkin-Huxley, but
% more like neurons in the cortex, being type-I.
% See Dayan and Abbott Sect 5.5, pp 166-172
% then Sect 6.1-2, pp.196-198 and Sect 6.6 p.224.
%
% For the original article see:
% Connor JA, Stevens CF, J Physiol 213:31-53 (1971)
%
% This adapted code is used to produce Figure 4.11 in the textbook:
% An Introductory Course in Computational Neuroscience,
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%% Parameters for plotting
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%% Basic cell properties for the model
V_L = -0.017;   % leak reversal potential
E_Na = 0.055;   % reversal for sodium channels
E_K = -0.072;   % reversal for potassium channels
E_A = -0.075;   % reversal for A-type current

G_L = 3e-8;     % specific leak conductance
G_Na = 1.2e-5;  % specific sodium conductance
G_K = 2e-6;     % specific potassium conductance
G_A = 4.77e-6;  % specific A-tpe potassium conductance

Cm = 0.1e-9;    % specific membrane capacitance

dt = 0.000005;  % time-step for the simulation

Vspike = 0.0;           % value for spike onset detection
Vendspike = -0.020;     % value for spike offset detection

% Part 1 for current step (A, C in figure), Part 2 for f-I curve
for figure_part = 1:2  
    
    if ( figure_part == 1 )     % current step
        tmin = -0.1;    % time of simulation commencement
        tmax=0.7;       % time to stop simulation
        
        Iappvec=0.85e-9;     % applied current is fixed
        
        istart = 0.1;   % time applied current starts
        ilength=0.5;    % length of applied current pulse
    else;               % f-I curve
        tmin = 0;       % time of simulation commencement
        tmax = 5;       % time to stop simulation
        istart = tmin;  % time applied current starts
        ilength = tmax-tmin;    % length of applied current pulse
        
        t=0:dt:tmax;    % time vector
        
        Iappvec = 0e-9:0.05e-9:3e-9; % set of applied current values
        
    end
    
    % nstart and nstop are indices of time vector 
    % for current onset and offset
    nstart = floor((istart-tmin)/dt)+1;
    nstop = floor((istart+ilength-tmin)/dt)+1;    
    Ibase = 0e-9;               % base current outside of step

    Ntrials = length(Iappvec);  % No. of points in f-I curve (1 in part 1)    
    rate = zeros(1,Ntrials);    % initialize for f-I curve
        
    t=tmin:dt:tmax;    % time vector
    
    for trial = 1:Ntrials;      % loop through trials (only once in part 1)
        
        Ie = Iappvec(trial)     % applied current value for this trial
        
        Iapp=Ibase*ones(size(t)); % Applied current, relevant in current-clamp mode
        for i= nstart:nstop % make non-zero for duration of current pulse
            Iapp(i) = Ie;
        end
        
        spike = zeros(size(t));     % vector to store spike times   
        inspike = 0;                % sets to 1 during a spike

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
        I_Na=zeros(size(t));    % to store sodium current vs time 
        I_K=zeros(size(t));     % to store potassium current vs time
        I_A=zeros(size(t));     % to store A-type current vs time
        I_L=zeros(size(t));     % to store leak current vs time
        
        for i = 2:length(t); % now see how things change through time
            
            Vm = V(i-1); % converts voltages to mV as needed in the equations on p.224 of Dayan/Abbott
            
            % find the steady state and time constant of all gating variables 
            % give the value of the membrane potential
            [m_inf, tau_m, h_inf, tau_h, n_inf, tau_n, a_inf, tau_a, ...
                b_inf, tau_b] = gating(Vm);
            
            m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
            h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
            n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
            
            a(i) = a(i-1) + (a_inf-a(i-1))*dt/tau_a;    % Update a
            b(i) = b(i-1) + (b_inf-b(i-1))*dt/tau_b;    % Update b
            
            I_L(i) = G_L*(V_L-V(i-1));                  % leak current
            
            I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i-1));   % sodium current
            
            I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i-1));      % potassium current
            
            I_A(i) = G_A*a(i)*a(i)*a(i)*b(i)*(E_A-V(i-1));      % A-type current
            
            Itot(i) = I_L(i)+I_Na(i)+I_K(i)+I_A(i)+Iapp(i);     % total current is sum of leak + active channels + applied current
            
            V(i) = V(i-1) + Itot(i)*dt/Cm;        % Update the membrane potential, V.
            
            % Section to record spike time on the upswing.
            % Note the use of the variable "inspike" that is set to zero
            % when the code is not in a spike allowing for spike detection.
            % Once the spike is detected, further time-points of high V are
            % not recorded as extra spikes until V drops and the variable
            % "inspike" is reset to zero.
            if ( inspike == 0 ) && ( V(i) > Vspike )    
                inspike = 1;                % Now "in a spike"
                spike(i) = 1;               % Record spike time
            else
                if ( inspike == 1 ) && ( V(i) < Vendspike )
                    inspike = 0;            % Allows detection of next spike
                end
            end
            
        end
        
        % Next line finds spike times in seconds from the array indices
        % produced by the function "find"
        spiketimes = dt*find(spike);
        if ( length(spiketimes) > 1 )   % if more than 1 spike
            isis = diff(spiketimes);    % find the inter-spike intervals
            rate(trial) = 1/isis(end);  % final rate is 1/ISI
        end
        
    end;        % Next trial for f-I curve
    
    %% Section to plot the appropriate curves for Fig 3.11
    if ( figure_part == 1 )             % single trial with current step
        figure(1)
        clf
        % Plot V vs time for current step
        subplot('Position',[0.12 0.6 0.36 0.36])
        plot(t,V,'k');
        xlabel('time, sec')
        ylabel('Membrane potential, V')
        axis([0 0.7 -0.075 0.06])
        
        % Zoom in on potassium currents around a spike
        subplot('Position',[0.12 0.12 0.36 0.36])
        plot(t,I_K*1e9,'k-');
        hold on
        plot(t,I_A*1e9,'k:')
        axis([0.2 0.24 -48 0])
        
        xlabel('Time (sec)')
        ylabel('Current (nA)')
        legend('I_{K}', 'I_{A}')
    else                                    % Many trials with different Iapp
        % Plot f-I curve
        subplot('Position',[0.6 0.6 0.36 0.36])
        plot(Iappvec*1e9,rate,'k');
        xlabel('Applied current (nA)')
        ylabel('Firing rate (Hz)')
    end
 
end;            % End loop through first two figure parts

%% Next show the V-dependence of the I_A gating variables
Vm = -0.100:0.001:0.040;        % Range of values for V
[m_inf, tau_m, h_inf, tau_h, n_inf, tau_n, a_inf, tau_a, ...
    b_inf, tau_b] = gating(Vm);

% Plot gating variables versus V
subplot('Position',[0.6 0.12 0.36 0.36])
plot(Vm*1000,a_inf,'k--');
xlabel('Membrane potential (mV)')
ylabel('Gating variable steady state')
hold on
plot(Vm*1000,b_inf,'k:');
plot(Vm*1000,n_inf,'k');
axis([-100 40 0 1])
legend('a_{\infty} ', 'b_{\infty} ', 'n_{\infty} ')

annotation('textbox',[0 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.5 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.0 0.52 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.5 0.52 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','D')

%% Notice the function is within the file containing the main code, but follows it here.
function[m_inf, tau_m, h_inf, tau_h, n_inf, tau_n, a_inf, tau_a, ...
    b_inf, tau_b] = gating(Vm);

% The function "gating" returns the steady states and the time constants of 
% the gating variables for the Connor-Stevens model.
% It should be sent the membrane potential, Vm, in units of Volts.
%
%
% Sodium and potassium gating variables are defined by the
% voltage-dependent transition rates between states, labeled alpha and
% beta. Written out from Dayan/Abbott, units are 1/ms.

alpha_m = 3.80e5*(Vm+0.0297)./(1-exp(-100*(Vm+0.0297)));
beta_m = 1.52e4*exp(-55.6*(Vm+0.0547));

alpha_h = 266*exp(-50*(Vm+0.048));
beta_h = 3800./(1+exp(-100*(Vm+0.018)));

alpha_n = 2e4*(Vm+0.0457)./(1-exp(-100*(Vm+0.0457)));
beta_n = 250*exp(-12.5*(Vm+0.0557));

% From the alpha and beta for each gating variable we find the steady
% state values (_inf) and the time constants (tau_) for each m,h and n.

tau_m = 1./(alpha_m+beta_m);      % time constant converted from ms to sec
m_inf = alpha_m./(alpha_m+beta_m);

tau_h = 1./(alpha_h+beta_h);      % time constant converted from ms to sec
h_inf = alpha_h./(alpha_h+beta_h);

tau_n = 1./(alpha_n+beta_n);      % time constant converted from ms to sec
n_inf = alpha_n./(alpha_n+beta_n);

% For the A-type current gating variables, instead of using alpha and
% beta, we just use the steady-state values a_inf and b_inf along with
% the time constants tau_a and tau_b that are found empirically
a_inf = (0.0761*exp(31.4*(Vm+0.09422))./(1+exp(34.6*(Vm+0.00117)))).^(1/3.0);
tau_a = 0.3632e-3 + 1.158e-3./(1+exp(49.7*(Vm+0.05596)));

b_inf = (1./(1+exp(68.8*(Vm+0.0533)))).^4;
tau_b = 1.24e-3 + 2.678e-3./(1+exp(62.4*(Vm+0.050)));

end