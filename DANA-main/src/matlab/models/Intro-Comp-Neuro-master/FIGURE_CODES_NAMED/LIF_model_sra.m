% LIF_model_sra.m
% Simulates a leaky integrate-and-fire neuron, which receives a single
% 200ms current pulse on each trial to cause spiking.
% Pulses of two amplitudes are used.
% The simulation includes an additional potassium conductance to provide 
% spike-rate adaptation. 
%
%
% This code was used to produce Figure 2.5 in the textbook:
% "An Introductory Course in Computational Neuroscience"
% by Paul Miller, Brandeis University 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;                      % clear all old variables from memory
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
Vth = -0.050;               % fixed threshold potential (to produce spike)
Vreset = -0.080;            % reset potential (post-spike)
Cm = 100e-12;               % total membrane capacitance
G_L = Cm/tau;               % total membrane conductance (leak conductance)

tsra = 0.200;               % time-constant for adaptation current
delta_G = 1e-9;             % potassium conductance increment per spike
Iapp = [240e-12 400e-12];   % Values of applied current to use
Ntrials = length(Iapp);     % Number of trials, one per applied current

%% Now loop through trials with different applied currents
for trial = 1:Ntrials;
    I = zeros(size(t));         % Initialize current to zero
    % Then generate the step current
    I(non:noff) = Iapp(trial)*ones(1,noff+1-non);
    V = E_L*ones(size(t));      % Initialize membrane potential
    spikes = zeros(size(t));    % Initialize vector of recorded spikes
    Gsra = zeros(size(t));      % Initialize adaptation conductance
    
    for i = 2:length(t)
        % Next line decays the adaptation conductance to zero        
        Gsra(i) = Gsra(i-1)*exp(-dt/tsra);
        % Next line calculates the steady state voltage directly from
        % the differential equation
        Vss = ( I(i) + G_L*E_L + Gsra(i)*E_K)/(G_L + Gsra(i));
        % Next line calculates the timescale for change in the membrane
        % potential
        taueff = Cm/(G_L+Gsra(i));
        % Next line updates the voltage using the exponential solution
        % for a first order linear ordinary differential equation
        V(i) = Vss + ( V(i-1)-Vss)*exp(-dt/taueff);
        
        if V(i) > Vth;              % if potential is above threshold
            spikes(i) = 1;          % record the spike at that time-point
            V(i) = Vreset;          % reset the potential
            Gsra(i) = Gsra(i) + delta_G;    % Update conductance for SRA
        end;
    end
    %% Integration through time is over
    % Now plot the membrane potential versus time
    
    figure(1)
    hold on;            % Prevents overwriting on the same subfigure
    % Next line positions each panel for each subfigure
    subplot('position',[0.45*trial-0.27 0.57 0.34 0.34])
    plot(t,1000*V,'k');      % Plot membrane potential versus time
    if ( trial == 1 )   % Label and mark y-axis for leftmost panels
        ylabel('V_m (mV)')
    end
    set(gca,'YTick',[-80 -50])
    axis([0 0.5 1000*(Vreset-0.005) 1000*(Vth+0.005)])
    if ( trial == 1)
        title('240pA')
    end
    if ( trial == 2)
        title('400pA')
    end
    
    % Lower subfigures plot adaptation conductance
    subplot('position',[0.44*trial-0.27 0.13 0.34 0.34])
    plot(t,Gsra*1e9,'k')
    xlabel('Time (sec)')
    if ( trial == 1 )
        ylabel('G_{SRA} (nS)')      % Y-label on left panel
    end
    
    axis([0 0.5 0 8e9*delta_G])     % Set x-range and y-range
    
end;        % End loop and go to next applied current

annotation('textbox',[0.00 0.97 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.00 0.53 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
