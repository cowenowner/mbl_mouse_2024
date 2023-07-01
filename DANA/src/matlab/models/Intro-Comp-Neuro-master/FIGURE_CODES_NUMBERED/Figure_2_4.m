% Figure_2_4.m
% Simulates a leaky integrate-and-fire neuron, which receives a single
% 200ms current pulse on each trial to cause spiking. 
% Pulses of two amplitudes are used for each of three models.
% Models differ in how the refractory conductance is included after a
% spike.
% 
% 
% This code was used to produce Figure 2.4 in the textbook:
% "An Introductory Course in Computational Neuroscience" 
% by Paul Miller, Brandeis University 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the default plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
    figure(1)
    clf;

%% Simulation parameters
dt = 0.0001;                % time-step
t = 0:dt:0.5;               % vector of time-points
ton = 0.15;                 % time to begin applied current (onset)
toff = 0.35;                % time to end applied current (offset)
non = round(ton/dt);        % time-point index of current onset
noff = round(toff/dt);      % time-point index of current offset

%% Parameters for the LIF neuron
tau = 0.010;                % membrane time constant
E_L = -0.070;               % leak potential (also resting potential)
E_K = -0.080;               % potassium reversal potential (for refractory conductance)
Vth0 = -0.050;              % baseline threshold potential (to produce spike)
Vreset = -0.080;            % reset potential (post-spike)
Cm = 100e-12;               % total membrane capacitance
G_L = Cm/tau;               % total membrane conductance (leak conductance)

tref = 0.002;               % Time constant of refractory period
delta_G = 500e-9;           % Increase in potassium conductance (method 2)
delta_Vth = 0.2;            % Increase in threshold voltage (method 3)

Iapp = [240e-12 400e-12];   % Two levels of applied current used per method
Ntrials = length(Iapp);     % Number of trials per method

%% Now loop through the three different methods for implementing a refractory period.
for tref_method = 1:3;
     
    % Loop through trials with different applied currents
    for trial = 1:Ntrials;      
        I = zeros(size(t));         % Initialize current to zero
        % Then generate the step current
        I(non:noff) = Iapp(trial)*ones(1,noff+1-non);
        V = E_L*ones(size(t));      % Initialize membrane potential
        spikes = zeros(size(t));    % Initialize vector of recorded spikes
        Gref = zeros(size(t));      % Initialize refractory conductance
        
        Vth = Vth0*ones(size(t));   % Initialize threshold voltage (which now can vary)
 
        % Ensure simulation does not commence in the refractory period of 
        % the last spike
        t_last_spike = -10*tref; 
        
        %% Commence integration through time
        for i = 2:length(t);
            % Next line decays the refractory conductance to zero
            Gref(i) = Gref(i-1)*exp(-dt/tref);
            % Next line decays the voltage threshold to its baseline value
            Vth(i) = Vth0 + (Vth(i-1)-Vth0)*exp(-dt/tref);
            % Next line calculates the steady state voltage directly from
            % the differential equation
            Vss = ( I(i) + G_L*E_L + Gref(i)*E_K)/(G_L + Gref(i));
            % Next line calculates the timescale for change in the membrane
            % potential
            taueff = Cm/(G_L+Gref(i));
            % Next line updates the voltage using the exponential solution 
            % for a first order linear ordinary differential equation
            V(i) = Vss + ( V(i-1)-Vss)*exp(-dt/taueff);
            
            % Now include the refractory period using one of three
            % different methods.
            if ( tref_method == 1 )     % Method (1) is forced voltage clamp
                if ( t(i) < t_last_spike + tref )   % If within tref of last spike
                    V(i) = Vreset;      % Clamp voltage to reset
                end
            end
            
            if V(i) > Vth(i)            % If membrane potential is above threshold
                spikes(i) = 1;          % Record a spike
                if (tref_method == 2 )  % Method (2) uses a potassium conductance
                        Gref(i) = Gref(i) + delta_G;    % Increase conductance
                else;                   % For methods (1) and (3) use reset
                    V(i) = Vreset;      % Reset the membrane potential
                    if (tref_method == 3)   % Method (3) raises threshold
                        Vth(i) = Vth0 + delta_Vth;  % Raise the threshold
                    end
                end
                t_last_spike = t(i);    % Keep track of time since spike
            end
        end;
        %% Integration through time is over
        % Now plot the membrane potential versus time
        
        figure(1)
        hold on;            % Prevents overwriting on the same subfigure
        % Next line positions each panel for each subfigure
        subplot('position',[0.45*trial-0.28 1.02-0.31*tref_method 0.34 0.22])
        plot(t,1000*V,'k');      % Plot membrane potential versus time
        hold on;            % Do not overwrite the subfigure
        plot(t,1000*Vth,'k--')   % Plot the threshold as a dashed line
        if ( trial == 1 )   % Label and mark y-axis for leftmost panels
            ylabel('V_m (mV)')
        end
        set(gca,'YTick',[-80 -40 0])
        axis([0.25 0.3 1000*(Vreset-0.005) 1000*(Vth0+0.05)]) % Zoom in on x-axis

        if ( tref_method == 1 )     % Give titles on top panels
            if ( trial == 1)
                title('240pA')
            end
            if ( trial == 2)
                title('400pA')
            end
        end
        set(gca,'XTick',[0.26 0.28])            % Set tick marks and
        set(gca,'XTickLabel',{'0.26', '0.28'})  % labels on bottom panel
        if ( tref_method == 3 )
            xlabel('Time (sec)')             
        end
        
    end
    
    if ( tref_method == 2 )     % Method with refractory conductance
        figure(2);              % Not included in book figure
        clf
        plot(t,Gref)            % Plot refractory conductance vs time
        figure(1)
    end
end

annotation('textbox',[0.00 0.97 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.00 0.66 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.00 0.35 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C')
