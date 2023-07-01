% Tutorial_8_4.m
%  Model of cerebellar conditioning from the output of an E-I chaotic
%  network to represent Granule Cells coupled to inhibitory Golgi cells.
%  The outputs of the Granule Cells connect to, via parallel fibers, and 
%  excite one Purkinje cell in this model.
%  The Purkinje cell also receives input from a Climbing Fiber, whose 
%  activity is treated only as a reinforcement signal to induce LTD 
%  in the active parallel fiber synapses. 
%  The output neuron of the model is a neuron of the Anterior Interpositus 
%  Nucleus, whose excitation depends on Granule Cell activity, whereas it 
%  is inhibited by the Purkinje Cell. Thus, a combination of Granule Cell
%  activity combined with a pause in Purkinje Cell activity produces
%  activity in the AIN cell, which leads to eyelid closure in this model.
%
%  The population of Granule Cells is activated by the Conditioned
%  Stimulus, which remains on through each trial. A later Unconditioned 
%  Stimulus causes Climbing Fiber activity, which is the reinforcement 
%  signal.
%
%  This code is a solution of Tutorial 8.4 of the textbook 
% "An Introductory Course in Computational Neuroscience" 
% by Paul Miller, Brandeis University 2016.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;                      % Clear all variables
dt = 0.0002;                % Time step for simulation
tmax = 2;                   % Maximum time of each trial
t = 0:dt:tmax;              % Vector of time points
Ntrials = 200;              % Number of trials for learning

NGcells = 200;              % No. of Granule Cells
NIcells = NGcells;          % No. of inhibitory (Golgi) cells

tCSon = 0.2;                % time of CS onset
lengthCS = tmax;            % length of CS
ICS0 = 2.0;                 % initial peak CS current

tUSon = 1.0;                % time for US onset
lengthUS = 0.020;           % length of US
IUS0 = 100.0;               % applied current for US

% Granule cell f-I curves:
IthGC = 10;                         % Thresholds of granule cells
IwidthGC = 1;                       % Width of f-I curve
rmaxGC = 100;                       % Max. firing rate
tauGC = 0.010;                      % membrane time constant
taus = 0.050;                       % base value for synaptic time constant 

% Inhibitory (Golgi) cell f-I curves
IthIC = 10;                 % fixed threshold
IwidthIC = 1;               % fixed width of sigmoid
rmaxIC = 100;               % fixed maximum rate 
tauIC = 0.005;              % fixed time constant

% Climbing fiber f-I curve
IthCF = 15;                 % fixed threshold
IwidthCF = 1;               % fixed width of sigmoid
rmaxCF = 30;               % fixed maximum rate 
tauCF = 0.002;              % fixed time constant

IthPC = 40;                 % fixed threshold
IwidthPC = 10;               % fixed width of sigmoid
rmaxPC = 30;                % fixed maximum rate 
tauPC = 0.010;              % fixed time constant

IthAIN = 25;                 % fixed threshold
IwidthAIN = 4;               % fixed width of sigmoid
rmaxAIN = 100;               % fixed maximum rate 
tauAIN = 0.010;              % fixed time constant

%% Now set all the connection strengths
%  Only WGP will vary via reinforcement learning

WGP0 = 10.0/NGcells;        % Base granule cell to Purkinje Cell connection strength
WGP = WGP0*ones(NGcells,1); % Vector of individual connections
dW0 = 0.1;                  % Scales rate of change of synaptic strength 

WGN = 5.0;                  % Connection from granule cell to nucleus (AIN)
WPN = 30.0;                 % Inhibitory strength from Purkinje Cell to nucleus (AIN)

% Produce the balanced chaotic network as in highD_chaos.m
% First the sets of random connection strengths
WGG = 2*rand(NGcells);          % Granule cells to granule cells connection
WGI = 2*rand(NGcells,NIcells);  % Granule cells to interneurons connection
WIG = 2*rand(NIcells,NGcells);  % Interneurons to granule cells connection
WII = 2*rand(NIcells);          % Interneuron to interneuron connections

% Next the sparse connection probabilities
connGG_prob = 0.05;             % Granule cells to granule cells connection
connIG_prob = 0.05;             % Granule cells to interneurons connection
connGI_prob = 0.05;             % Interneurons to granule cells connection
connII_prob = 0.05;             % Interneuron to interneuron connections

% Binary connectivitiy matrices produced based on probability of existence
% of a connection
connGG = rand(NGcells)<connGG_prob;
connIG = rand(NIcells,NGcells)<connIG_prob;
connGI = rand(NGcells,NIcells)<connGI_prob;
connII = rand(NIcells)<connII_prob;

% Full connection weight matrix, containing the random numbers but with
% sparse connections
WGG = WGG.*connGG;
WIG = WIG.*connIG;
WGI = WGI.*connGI;
WII = WII.*connII;

% Ensure balance at the individual cell level so that each neuron receives
% the same total excitatory connection strength and same total inhibitory
% connection strength as all others.
for i = 1:NGcells
    WGG(:,i) = WGG(:,i)/mean(WGG(:,i));
    WIG(:,i) = WIG(:,i)/mean(WIG(:,i));
end
for i = 1:NIcells
    WGI(:,i) = WGI(:,i)/mean(WGI(:,i));
    WII(:,i) = WII(:,i)/mean(WII(:,i));
end
WGG0 = 90/NGcells;          % Mean magnitude of connection strengths E-to-E
WGI0 = 80/NGcells;          % Mean magnitude of connection strengths E-to-I
WIG0 = 90/NIcells;          % Mean magnitude of connection strengths I-to-E
WII0 = 80/NIcells;          % Mean magnitude of connection strengths I-to-I

WGG = WGG*WGG0;             % Final connection strength matrix E-to-E
WGI = WGI*WGI0;             % Final connection strength matrix E-to-I
WIG = WIG*WIG0;             % Final connection strength matrix I-to-E
WII = WII*WII0;             % Final connection strength matrix I-to-I

%% Now set the vectors of currents for the conditioned stimulus (ICS) and 
%   the unconditioned stimulus (IUS)
ICS = zeros(size(t));       % set up input current vectors, one for CS
IUS = zeros(size(t));           % and one for US
for i = floor(tCSon/dt)+1:floor((tCSon+lengthCS)/dt)+1
    ICS(i) = ICS0;              % Conditioned Stimulus is on
end
for i = floor(tUSon/dt)+1:floor((tUSon+lengthUS)/dt)+1
    IUS(i) = IUS0;              % Unconditioned Stimulus is on
end

%% Loop through all of the trials
for trial = 1:Ntrials; trial
    rGC = zeros(length(t),NGcells);     % rates of all Granule Cells
    S = zeros(length(t),NGcells);       % synaptic output of granule cells
    rI = zeros(length(t),NIcells);      % rates of all inhibitory Golgi cells
    rPC = zeros(length(t),1);           % rate of Purkinje cell
    rCF = zeros(length(t),1);           % rate of spikes in Climbing Fiber
    rAIN = zeros(length(t),1);          % firing rate of Interpositus Nucleus
    ItotPC = zeros(size(t));            % current to Purkinje cells

    
    for i = 2:length(t)                 % Loop through time in one trial
        % Evaluate inputs and rates of granule cells        
        ItotGC = rGC(i-1,:)*WGG - rI(i-1,:)*WIG + ICS(i-1);
        rinfGC = rmaxGC./(1+exp(-(ItotGC-IthGC)/IwidthGC));    % rinf is a vector of the steady state firing rate                           
        rGC(i,:) = rinfGC + (rGC(i-1,:)-rinfGC)*exp(-dt/tauGC);   % update r from r(i-1) for all cells

        % Total synaptic input from Granule Cells to Purkinje Cell 
        % will be proportional to firing rate of Granule Cell
        Sinf = rGC(i,:)/rmaxGC;
        S(i,:) = Sinf + (S(i-1,:)-Sinf).*exp(-dt./taus);

        % Now evaluate inputs and rates of inhibitory (Golgi) cells
        ItotIC = rGC(i-1,:)*WGI- rI(i-1,:)*WII;       
        rinfI = rmaxIC./(1+exp(-(ItotIC-IthIC)/IwidthIC));    % rinf is a vector of the steady state firing rate
        rI(i,:) = rinfI + (rI(i-1,:)-rinfI)*exp(-dt/tauIC);   % update r from r(i-1) for all cells
        
        % Now evaluate inputs and rate of climbing fiber
        ItotCF = IUS(i);
        rinfCF = rmaxCF./(1+exp(-(ItotCF-IthCF)/IwidthCF)).^2;
        rCF(i) = rinfCF + (rCF(i-1)-rinfCF)*exp(-dt/tauCF);
        
        % Now evaluate inputs and rate of the Purkinje cell
        ItotPC(i) = rGC(i-1,:)*WGP;
        rinfPC = rmaxPC./(1+exp(-(ItotPC(i)-IthPC)/IwidthPC)).^2;
        rPC(i) = rinfPC + (rPC(i-1)-rinfPC)*exp(-dt/tauPC);

        % Now evaluate inputs and rate of the interpositus nuceleus
        ItotAIN = WGN*mean(rGC(i-1,:)) - rPC(i-1)*WPN; % excitation by mossy fibers and inhibition by Purkinje Cell
        rinfNI = rmaxAIN./(1+exp(-(ItotAIN-IthAIN)/IwidthAIN)).^2;
        rAIN(i) = rinfNI + (rAIN(i-1)-rinfNI)*exp(-dt/tauAIN);
    end   % end time loop
    
    % Now to adjust the synapses from GC to PC 
    ratecorr = zeros(NGcells,1);
    for cell = 1:NGcells
        % coincidence between granule cell output and climbing fiber spikes
        ratecorr(cell) = sum(S(:,cell).*rCF)*dt;     
    end
   
    ms = mean(S,1);                     % Mean synaptic output for each granule cell
    dW = -dW0*(ratecorr-0.1);   % LTD for correlated firing, otherwise LTP
    WGP = WGP.*(1.0+dW);                % Update all granule to PC synapses
    WGP = max(WGP,0.0);                 % Strength can not be less than zero
    WGP = min(WGP,5*WGP0);              % Strength has a maximum value

    % In the next line record the minimum input to Purkinje Cells in the
    % period when the conditioned response should arise
    min_IPC(trial) = min(ItotPC(floor((tCSon+0.1)/dt):end));
    
    %% Now plot the figures on the first trial and later trials
    % Plot the chaotic activity in trial 1 (remains unchanged)
    if ( trial == 1 ); figure(1); subplot(2,2,1); plot(t,rGC); axis([ 0 tmax -5 105 ]); end
    
    % Every 25 trials update the figure
    if ( ( trial == 1 ) || ( mod(trial,25) == 0 ) )
        min_input = min_IPC(trial)      % Output this information
        % Plot the AIN rate, which is high for a conditioned response
        figure(1); subplot(2,2,2); cla; plot(t,rAIN);  hold on;

        % Plot the Purkinje Cell rate, which should dip at the conditioned
        % response.
        plot(t,rPC,'r'); xlabel('time'); ylabel('Firing rate (red=Purkinje; blue=Nucleus)');
        
        % Plot the weights with plasticity (Granule Cell to Purkinje)
        figure(1); subplot(2,2,3); cla; plot(WGP); axis([0 NGcells 0 5*WGP0 ])
        xlabel('Granule Cell Index'); ylabel('Connection Strength, WGP'); hold on;

        % Plot the current input to the Purkinje Cell from Granule Cells
        figure(1); subplot(2,2,4); cla; plot(t,ItotPC); axis([0 tmax 0 1200])
        xlabel('time'); ylabel('Purkinje Cell Input');
        drawnow
    end
end  % end of trial loop