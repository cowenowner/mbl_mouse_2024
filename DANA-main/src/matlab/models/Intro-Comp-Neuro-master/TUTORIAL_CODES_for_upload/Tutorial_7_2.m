% Tutorial_7_2.m
%
% This code produces those dynamical regimes possible with two units,
% with a firing-rate model of a neural circuit.
% Two units (representing neural
% assemblies) are simulated with firing rates that increase linearly with
% input current above a system-dependent threshold (Ithresh) and that saturate at a value
% rmax = 100Hz.
% Units can excite or inhibit each other to different degrees according to
% the connectivity matrix W in which W(i,j) represents connection strength
% from unit j to unit i.
% The flag "attractor-flag" is used to determine the type of simulation.
%
% This code is a solution of Tutorial 7.2 in the textbook,
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University, 2017
%%

clear

question_number = 9;         % Allowed values are 1, 2, 3, 5, 6
% question_number determines the connection strengths and the type of
% dynamical system produced.
% 1: point attractor                            Figure 1 in the paper
% 2: multistable point attractor system         Figure 2 in the paper
% 3: inhibition stabilized point attractor      Figure 3 in the paper
% 4: should not be used with two units
% 5: marginal state (line attractor)            Figure 5 in the paper
% 6: oscillator (bistable)                      Figure 6 in the paper
% 7: should not be used with two units
% 8: should not be used with two units
% 9: should not be used with two units
tmax = 3;  % default value of maximum time to simulate
Ntrials = 1;

switch question_number
    case 1,
        % case 1 is a point attractor system. Upon input the activity of
        % the circuit changes but to a new, single point attractor -- rates
        % are determined by the input and not history-dependent and become
        % time-independent. Any noise would merely produce fluctuations
        % about the stable state.
        Ithresh = [-5; -10];
        W = [0.6 1 ;-0.2 0];
        rinit1 = [10; 50];
        Iapp1 = [10; 0];
        Iapp2 = [30; 0];
        Nunits = 2;             % Number of units to simulate
        
        
    case 2,
        % case 2 is a multistable system with three point attractor states
        % (all cells at low rates or cell 1 high and cell 3 low or cell 3
        % high and cell 1 low. In such a system some cells can produce more
        % graded activity (cell 2) due to changes in input from the others.
        % Initial conditions or the types of input pulse determine the
        % final attractor state, so the system exhibits memory.
        Ithresh = [10; 5];
        W = [1.2 -0.3; -0.2 1.1];
        rinit1 = [50; 55];
        Iapp1 = [ 0; 50];
        Iapp2 = [50; 0];
        Nunits = 2;             % Number of units to simulate
        
        
    case 3,
        % case 3 is a point attractor in the inhibition stabilized regime
        % such that application of a negative current pulse to the
        % inhibitory cell 2 produces an eventual increase in firing rate
        % of all cells, including the one receiving external inhibition.
        Ithresh = [-10; 0];
        W = [2.5 2; -3.0 -2];
        rinit1 = [40; 0];
        Iapp1 = [0; -5];
        Iapp2 = [0; -10];
        Nunits = 2;             % Number of units to simulate
        
        
    case 4,
        % case 4 produces a marginal state, or line attractor as the line
        % described by r1 + r2 = 200/3 is a continuous set of fixed points
        % in the range 0 < r1,r2 < rmax (= 100).
        % Final set of activities depends on the initial values and current
        % pulses are integrated.
        Ithresh = [-10; -10];
        W = [0.8 -0.2; -0.4 0.6];
        rinit1 = [30; 75];
        Iapp1 = [1; 0];
        Iapp2 = [2; 0];
        Nunits = 2;             % Number of units to simulate
        
        
    case 5,
        % case 6 produces an oscillator with a triphasic rhythm. Changes in
        % initial conditions alter the phase of the oscillation but not the
        % pattern, which is an orbit attractor or limit cycle.
        Ithresh = [0; 20];
        W = [2.0 1; -1.5 0];
        rinit1 = [80; 0];
        Iapp1 = [0; 0];
        Iapp2 = [-10; 0];
        Nunits = 2;             % Number of units to simulate
        
    case 6,
        % case 4 has two point attractors in the inhibition stabilized
        % regime. The system is bistable and can be switched between
        % states of low firing rate.
        Ithresh = [-10; -5; 5];
        W = [1.5 0 1; 0 2 1 ; -2.5 -3 -1 ];
        rinit1 = [10; 5; 10];
        rinit2 = [30; 5; 5];
        Iapp1 = [20; 0; 0];
        Iapp2 = [0; 25; 0];
        Nunits = 3;
        
    case 7,
        % case 7 produces two distinct oscillators of different
        % frequencies. Changes in initial conditions -- or a current pulse
        % -- can switch between oscillators.
        Ithresh = [-18; -15; 0];
        W = [2.2 -0.5 0.9; -0.7 2 1.2; -1.6 -1.2 0];
        rinit1 = [50; 15; 0];
        rinit2 = [10; 15; 20];
        Iapp1 = [0; 60; 0];
        Iapp2 = [50; 0; 0];
        Imax = 55;
        Ntrials = 1;
        Nunits = 3;
        
    case 8,
        %  case 8 is chaotic. A small change in initial conditions causes
        %  changes in the behavior that accumulate over time, so that for
        %  example the timings of small peaks of activity in r1 eventually
        %  become unpredictable. Note the near periodicity of chaotic
        %  activity that can produce peaks in the power spectrum.
        Ithresh = [-10; -20; 10 ];
        W = [2.05 -0.2  1.2  ; -0.05  2.1   0.5 ; -1.6 -4  0 ];
        rinit1 = [85; 15; 10];
        rinit2 = [86; 15; 10];
        Iapp1 = [0; 0; 0];
        Iapp2 = [0; 0; 0];
        tmax = 3;
        Ntrials = 2;
        Nunits = 3;
        
    case 9,
        %  case 9 is a heteroclinic with saddle points at (100,0,0) and (0,100,100).
        %  orbits move to the vicinity of these saddle points, so timescale
        %  of initial transient depends on initial conditions and how
        %  closely the saddle points are approached.
        Ithresh = [-2; -1; -1];
        W = [0.98  -0.015 -0.01; 0 0.99 -0.02 ; -0.02 0.005 1.01];
        rinit1 = [0; 95; 99.9];
        rinit2 = [0; 95; 99.8];
        Iapp1 = [0; 0; 0];
        Iapp2 = [0; 0; 0];
        tmax = 20;
        Ntrials = 2;
        Nunits = 3;
        
    otherwise
        disp('system undefined')
end

%% Figure positions for subplots
subfigA = [0.14 0.79 0.78 0.16];
subfigB = [0.14 0.45 0.78 0.16];
subfigC = [0.14 0.11 0.78 0.16];
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');


%% Set up the time vector
dt = 0.0001;
tvec = 0:dt:tmax;
Nt = length(tvec);

r = zeros(Nunits,Nt);   % array of rate of each cell as a function of time
rmax = 100;             % maximum firing rate
tau = 0.010;            % base time constant for changes of rate

%% Set up details of an applied current used in some systems
Iapp = zeros(Nunits,Nt);    % Array of time-dependent and unit-dependent current
Ion1 = 1;                    % Time to switch on
Ion2 = 2;                    % Time to switch on
Idur = 0.2;                 % Duration of current

non1 = round(Ion1/dt);            % Time step to switch on current
noff1 = round((Ion1+Idur)/dt);    % Time step to switch off current
non2 = round(Ion2/dt);            % Time step to switch on current
noff2 = round((Ion2+Idur)/dt);    % Time step to switch off current

% Add the applied current pulses
Iapp(:,non1:noff1) = Iapp1*ones(1,noff1-non1+1);
Iapp(:,non2:noff2) = Iapp2*ones(1,noff2-non2+1);

%% Set up axes and label them for the figures

figure(question_number)          % Separate figure for each circuit
clf
subplot('position',subfigA)     % Subplot for Applied current
axis([0 tmax -15 35])
hold on
xlabel('Time, sec')
ylabel('Applied Current')

subplot('position',subfigB)     % Subplot for firing rate response
axis([0 tmax -5 rmax+5])
hold on
xlabel('Time, sec')
ylabel('Firing rate, Hz')

for trial = 1:Ntrials
    if ( trial == 1 )
        r(:,1) = rinit1;                % Initialize firing rate
    else
        r(:,1) = rinit2;                % Initialize firing rate
    end
    %% Now simulate through time
    for i=2:Nt
        I = W'*r(:,i-1) + Iapp(:,i-1);                   % total current to each unit
        newr = r(:,i-1) + dt/tau*(I-Ithresh-r(:,i-1));  % Euler-Mayamara update of rates
        r(:,i) = max(newr,0);                           % rates are not negative
        r(:,i) = min(r(:,i),rmax);                      % rates can not be above rmax
        
        %% Now plot figures rate vs time for each color-coded unit in figure 1,
        %  with a new panel for each trial updating every 1/10 of the
        %  simulation
        
        if ( mod(i,round(Nt/20) ) == 0 )
            
            subplot('position',subfigA)
            plot(tvec(1:i),Iapp(1,1:i),'k-')
            hold on
            plot(tvec(1:i),Iapp(2,1:i),'r-')
            
            if ( trial == 1 )
                subplot('position',subfigB)
            else
                subplot('position',subfigC)
            end
            plot(tvec(1:i),r(1,1:i),'k')
            hold on
            plot(tvec(1:i),r(2,1:i),'r')
            
            pause(0.05); drawnow
        end
    end
    
end
