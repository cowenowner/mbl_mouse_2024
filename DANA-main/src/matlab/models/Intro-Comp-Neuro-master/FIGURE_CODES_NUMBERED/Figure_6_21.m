% Figure_6_21.m 
% Code for producing a model of the head direction 
% system,which operates by integration of an angular-velocity input.
%
%
%   The bump moves because the inhibitory feedback is directional, with two
%   competing inhibitory populations biasing feedback in one direction or 
%   the other. Asymmetric input enhances one of the population's activity 
%   to generate drift.
%
%   Based on papers by Zhang (1998) and Song & Wang (2005).
%
%   This code is used to produce Figure 6.21 in Chapter 6 of the book
%
%   "An Introductory Course in Computational Neuroscience"
%
%   by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;                          % Clear all variables

question_part = 'A';            % Can be 'A' or 'B' for "EI + IE" or "II:
stim_random_on = 0;             % Set to 1 to include randomness in the stimulus
dt=0.0005;                      % Time step in sec (0.1ms)
tmax = 20;                      % Maximum time for integration 
t=0:dt:tmax;                    % Create time vector

Ncells = 50;                    % Number of cells around the ring

rE=zeros(length(t),Ncells);     % Rate matrix for all timesteps and all cells
rI1=zeros(length(t),Ncells);     % Rate matrix for all timesteps and all cells
rI2=zeros(length(t),Ncells);     % Rate matrix for all timesteps and all cells
hE=zeros(1,Ncells);             % Applied thalamic input to all E cells
hI=zeros(1,Ncells);             % Applied thalamic input to all I cells
IstimE=zeros(1,Ncells);         % Total current to all E cells
IstimI=zeros(1,Ncells);         % Total current to all I cells

f_of_I = @(x) (x>0).*x;

%% Cortical network parameters
switch question_part
    case 'A'                % case 'A' is akin to the standard EI ring network
        I0E = 4;            % Background input to E cells minus threshold
        I0I = -25;          % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 8;           % Mean E to E connection weight
        WEI0 = 6;           % Mean E to I connection weight
        WIE0 = -2;          % Mean I to E connection weight
        WIE1shift = pi/5;       % Difference in tuning preference for I cells connecting to E cells
        WIE2shift = -pi/5;      % Difference in tuning preference for I cells connecting to E cells
        WII0 = 0;             % Mean I to E connection weight
        WII11shift = 0;       % Difference in tuning preference for I cells connecting to I cells
        WII21shift = 0;       % Difference in tuning preference for I cells connecting to I cells
        WII12shift = 0;       % Difference in tuning preference for I cells connecting to I cells
        WII22shift = 0;       % Difference in tuning preference for I cells connecting to I cells
    
    case 'B'                % case 'B' is a purely inhibitory network
        I0E = -1;           % Background input to E cells minus threshold
        I0I = 25;           % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 0;           % Mean E to E connection weight (no E-Cells)
        WEI0 = 0;           % Mean E to I connection weight (no E-Cells)
        WIE0 = 0;           % Mean I to E connection weight (no E-Cells)
        WIE1shift = -pi/10;       % Difference in tuning preference for I cells connecting to E cells
        WIE2shift = pi/10;        % Difference in tuning preference for I cells connecting to E cells
        WII0 = -20;             % Mean I to E connection weight
        WII11shift = pi-pi/20;  % Difference in tuning preference for I cells connecting to I cells
        WII12shift = pi;        % Difference in tuning preference for I cells connecting to I cells
        WII21shift = pi;        % Difference in tuning preference for I cells connecting to I cells
        WII22shift = pi+pi/20;  % Difference in tuning preference for I cells connecting to I cells
end

% Now produce all the within-network connections as cosine functions. Note
% that 1 is added to the cosine so the variation of the term within
% parentheses has a minimum of 0 and a maximum of 2.
% The initial terms WEE0, WEI0, and WIE0 correspond then to the mean weight
% of that type of connection.
% The I-to-E connection has an extra "WIEshift" term which can be set to pi
% so that inhibition is from cells with opposite LGN input.
for cell1 = 1:Ncells
    for cell2 = 1:Ncells
        WEE(cell1,cell2) = WEE0*(1+cos(2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WEI(cell1,cell2) = WEI0*(1+cos(2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WIE1(cell1,cell2) = 0.5*WIE0*(1+cos(WIE1shift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WIE2(cell1,cell2) = 0.5*WIE0*(1+cos(WIE2shift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WII11(cell1,cell2) = 0.5*WII0*(1+cos(WII11shift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WII12(cell1,cell2) = 0.5*WII0*(1+cos(WII12shift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WII21(cell1,cell2) = 0.5*WII0*(1+cos(WII21shift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WII22(cell1,cell2) = 0.5*WII0*(1+cos(WII22shift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
    end
end

%% Stimulus input parameters
epsilon = 1;              % Variation of input across cells from 1-epsilon to 1+epsilon
cuecell = Ncells/2;         % Cell with peak of the LGN input
hE = zeros(1,Ncells);       % Vector for inputs to each E cell when input is on
hI = zeros(1,Ncells);       % Vector for inputs to each I cell when input is on
hstart = tmax/4;               % Time to begin input
hend = tmax/2;                % Time to end input
velocities = [0:40];        % Range of inputs to use, representing angular velocity
Ntrials = length(velocities);  % Number of trials
noise_amp = 0.2;            % Noise heterogeneity in the input when stim_random_on=1
trajectory = zeros(length(t),Ntrials);
%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%% Now loop through the set of different input biases representing velocity
for trial = 1:Ntrials
    velocity = velocities(trial);       % Input velocity to be used
    
    rE(1,cuecell) = 50;
    rI1(1,cuecell) = 20;
    rI2(1,cuecell) = 20;
    % Begin the time integration
    for i = 2:length(t)
        
        IstimE = 0;                     % External input to E-cells
        IstimI1 = velocity;             % External input to I-cells (1)
        IstimI2 = - velocity;           % External input to I-cells (2)
        
        % Update rates of all excitatory units based on rates of all units
        % at the previous time point
        rE(i,:) = rE(i-1,:)*(1-dt/tauE) + ...
            dt*f_of_I(IstimE + rE(i-1,:)*WEE + rI1(i-1,:)*WIE1 + rI2(i-1,:)*WIE2 + I0E)/tauE;
        
        % Update rates of all inhibitory units based on rates of all units
        % at the previous time point
        rI1(i,:) = rI1(i-1,:)*(1-dt/tauI) + ...
            dt*f_of_I(IstimI1 + rE(i-1,:)*WEI + ...
            rI1(i-1,:)*WII11 +rI2(i-1,:)*WII21 + I0I)/tauI;
        rI2(i,:) = rI2(i-1,:)*(1-dt/tauI) + ...
            dt*f_of_I(IstimI2 + rE(i-1,:)*WEI + ...
            rI1(i-1,:)*WII12 +rI2(i-1,:)*WII22+I0I)/tauI;
        
        rE(i,:) = max(rE(i,:),0);            % Rates cannot be less than 0
        rI1(i,:) = max(rI1(i,:),0);          % Rates cannot be less than 0
        rI2(i,:) = max(rI2(i,:),0);          % Rates cannot be less than 0
        
    end
    
    %% Now plot the results for the contrast used
    figure(1)
    if ( trial == 1 )
        clf;                    % Clear figure on first trial
    end
    
    subplot(2,1,1)
    % plot rate of all E cells at end of simulation
    plot(rE(length(t),:),'g')
    hold on
    xlabel('cell index')
    ylabel('rate of E-neurons')
    
    subplot(2,1,2)
    % plot rate of all I cells at end of simulation
    plot(rI1(length(t),:),'g')
    hold on
    plot(rI2(length(t),:),'r')
    xlabel('cell index')
    ylabel('rate of I-neurons')
    
    figure(2)
    if ( trial == 1 )
        clf;                    % Clear figure on first trial
    end
    subplot(2,1,1)
    % plot rate of excitatory cuecell as a function of time
    plot(t,rE(:,cuecell))
    hold on
    % plot rate of E cell with null direction input as a function of time
    plot(t,rE(:,1+mod(cuecell+Ncells/2-1,Ncells)))
    xlabel('time')
    ylabel('rate of excitatory cell')
    
    
    subplot(2,1,2)
    % plot rate of inhibitory cuecell as a function of time
    plot(t,rI1(:,cuecell))
    hold on
    plot(t,rI2(:,cuecell))
    % plot rate of I cell with null direction input as a function of time
    plot(t,rI1(:,1+mod(cuecell+Ncells/2-1,Ncells)))
    xlabel('time')
    ylabel('rate of inhibitory cell')
    
    % Now extract the peak location of the bump as a function of time using 
    % the max function and store it in the array "traj"
    switch question_part
        case 'A'
            [vals traj(:,trial) ] = max(rE,[],2);
        case 'B'
            [vals traj(:,trial) ] = max(rI1,[],2);
    end
    figure(21)
    plot(t,traj(:,trial))
    hold on
    drawnow
end % Next velocity

v = diff(traj);         % measured velocities based in difference in peak locations

% The following line only works for positive velocities. 
% If velocities are negative the line should read:
% v = -mod(-v,Ncells);
% The mod sign is taken to account for a shift of the index of the bump
% center from the last to the first (or the first to the last)
% corresponding to a shift of +1 (or -1).
v = mod(v,Ncells);                  % Because angular position repeats  

meanv = mean(v(end/2:end,:))/dt/Ncells; % The mean velocity in cycles per second

% Now plot measured "output velocity as a function of "input velocity"
figure(22)
switch question_part
    case 'A'
        plot(velocities,meanv,'ok')
        hold on
    case 'B'
        plot(velocities,meanv,'*k')
        hold on
end

xlabel('Input Bias')
ylabel('Angular Velocity (rev/s)')
axis([0 40 0 6])