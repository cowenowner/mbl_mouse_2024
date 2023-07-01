% Figure_6_20.m
% 
%  Firing-rate model code for producing a stable "bump" of activity on a
%  ring, based on a small variation of the orientation selectivity network.
%
%  The bump can move if the inhibitory feedback is directional, in which
%  case it can represent head direction cells.
%  See paper by Compte et al Cerebral Cortex (2000) for a spiking neuron
%  model.
%
%  This code was used to produce Figure 6.19 in the book:
%  An Introductory Course in Computational Neuroscience
%  by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;                          % Clear all variables

question_part = 'C';            % Can be 'A' or 'B' or 'C' or 'D' or 'E'
stim_random_on = 0;             % Set to 1 to include randomness in the stimulus
dt=0.0001;                      % Time step in sec (0.1ms)
tmax =8;                      % Time step in sec (300ms)
t=0:dt:tmax;                    % Create time vector
Nt = length(t);

Ncells = 50;                    % Number of cells around the ring

rE=zeros(Nt,Ncells);     % Rate matrix for all timesteps and all cells
rI=zeros(Nt,Ncells);     % Rate matrix for all timesteps and all cells
hE=zeros(1,Ncells);             % Applied thalamic input to all E cells
hI=zeros(1,Ncells);             % Applied thalamic input to all I cells
IstimE=zeros(1,Ncells);         % Total current to all E cells
IstimI=zeros(1,Ncells);         % Total current to all I cells

%% Cortical network parameters
switch question_part
    % A: % Stationary bump based on excitatory feedback. Activity can be on or off.
    case 'A'
        AE = 2;             % Maximum LGN input to E cells
        AI = 0;             % Maximum LGN input to I cells
        I0E = -2;           % Background input to E cells minus threshold
        I0I = -4;           % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 8;           % Mean E to E connection weight
        WEI0 = 3;           % Mean E to I connection weight
        WIE0 = -3;          % Mean I to E connection weight
        WIEshift = 0.0;       % Difference in tuning preference for I cells connecting to E cells
        WII0 = 0;           % Mean I to E connection weight
        WIIshift = 0;       % Difference in tuning preference for I cells connecting to E cells
        
        % B: Stationary bump based on excitatory feedback. Activity is always on.
    case 'B'
        AE = 10;            % Maximum LGN input to E cells
        AI = 0;             % Maximum LGN input to I cells
        I0E = 1;            % Background input to E cells minus threshold
        I0I = -4;           % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 8;           % Mean E to E connection weight
        WEI0 = 3;           % Mean E to I connection weight
        WIE0 = -3;          % Mean I to E connection weight
        WIEshift = 0;       % Difference in tuning preference for I cells connecting to E cells
        WII0 = 0;           % Mean I to E connection weight
        WIIshift = 0;       % Difference in tuning preference for I cells connecting to E cells
    case 'B2'
        AE = 10;            % Maximum LGN input to E cells
        AI = 0;             % Maximum LGN input to I cells
        I0E = 2;            % Background input to E cells minus threshold
        I0I = -1;           % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 6;           % Mean E to E connection weight
        WEI0 = 3;           % Mean E to I connection weight
        WIE0 = -2;          % Mean I to E connection weight
        WIEshift = 0;       % Difference in tuning preference for I cells connecting to E cells
        WII0 = 0;           % Mean I to E connection weight
        WIIshift = 0;       % Difference in tuning preference for I cells connecting to E cells
        
        % C: Moving bump, symmetric excitatory and asymmetric inhibitory feedback.
    case 'C'
        AE = 0;            % Maximum LGN input to E cells
        AI = 0;             % Maximum LGN input to I cells
        I0E = 2;            % Background input to E cells minus threshold
        I0I = -1;           % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 6;           % Mean E to E connection weight
        WEI0 = 3;           % Mean E to I connection weight
        WIE0 = -2;          % Mean I to E connection weight
        WIEshift = 0.05*pi;       % Difference in tuning preference for I cells connecting to E cells
        WII0 = 0;           % Mean I to E connection weight
        WIIshift = 0;       % Difference in tuning preference for I cells connecting to E cells
    case 'C2'
        AE = 0;            % Maximum LGN input to E cells
        AI = 0;             % Maximum LGN input to I cells
        I0E = 4;            % Background input to E cells minus threshold
        I0I = -1;           % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 8;           % Mean E to E connection weight
        WEI0 = 3;           % Mean E to I connection weight
        WIE0 = -3;          % Mean I to E connection weight
        WIEshift = 0.05*pi;       % Difference in tuning preference for I cells connecting to E cells
        WII0 = 0;           % Mean I to E connection weight
        WIIshift = 0;       % Difference in tuning preference for I cells connecting to E cells
        
        % D: Stationary bump attractor based on symmetric I-to-I feedback.
    case 'D'
        AE = 0;             % Maximum LGN input to E cells
        AI = 20;            % Maximum LGN input to I cells
        I0E = -10;          % Background input to E cells minus threshold
        I0I = 20;           % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 0;           % Mean E to E connection weight
        WEI0 = 0;           % Mean E to I connection weight
        WIE0 = 0;           % Mean I to E connection weight
        WIEshift = 0;       % Difference in tuning preference for I cells connecting to E cells
        WII0 = -10;           % Mean I to E connection weight
        WIIshift = pi;       % Difference in tuning preference for I cells connecting to E cells
        
        % E: A moving bump attractor based on asymmetric I-to-I feedback.
    case 'E'
        AE = 0;            % Maximum LGN input to E cells
        AI = 20;             % Maximum LGN input to I cells
        I0E = -10;            % Background input to E cells minus threshold
        I0I = 20;           % Background input to I cells minus threshold
        tauE = 0.050;       % Time constant for E cells
        tauI = 0.005;       % Time constant for I cells
        WEE0 = 0;           % Mean E to E connection weight
        WEI0 = 0;           % Mean E to I connection weight
        WIE0 = 0;          % Mean I to E connection weight
        WIEshift = 0.05*pi;       % Difference in tuning preference for I cells connecting to E cells
        WII0 = -10;           % Mean I to E connection weight
        WIIshift = pi+0.001*pi;       % Difference in tuning preference for I cells connecting to E cells
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
        WIE(cell1,cell2) = WIE0*(1+cos(WIEshift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
        WII(cell1,cell2) = WII0*(1+cos(WIIshift+2*pi*(cell2-cell1)/Ncells) ) / Ncells;
    end
end

%% Stimulus input parameters
epsilon = 0.1;              % Variation of input across cells from 1-epsilon to 1+epsilon
cuecell = Ncells/2;         % Cell with peak of the LGN input
hE = zeros(1,Ncells);       % Vector for inputs to each E cell when input is on
hI = zeros(1,Ncells);       % Vector for inputs to each I cell when input is on
hstart = tmax/4;               % Time to begin input
hend = tmax/2;                % Time to end input
contrasts = [1];    % Range of contrasts to use
noise_amp = 0.2;            % Noise heterogeneity in the input when stim_random_on=1
%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

% Now loop through the set of different contrasts
for trial = 1:length(contrasts)
    c = contrasts(trial);               % Contrast to be used
    IstimE = zeros(1,Ncells);           % Define vector for input current to E cells
    IstimI = zeros(1,Ncells);           % Define vector for input current to E cells
    
    % Now set the input current that varies across cells
    for cell = 1:Ncells
        hE(cell) = AE*c*(1 + epsilon*cos(2*pi*(cell-cuecell)/Ncells));
        hI(cell) = AI*c*(1 + epsilon*cos(2*pi*(cell-cuecell)/Ncells));
    end
    if ( stim_random_on )
        hE = hE.*(1+noise_amp*(rand(1,Ncells)-0.5));
        hI = hI.*(1+noise_amp*(rand(1,Ncells)-0.5));
    end
    rE(1,10) = 0.1;
    % Begin the time integration
    for i = 2:Nt
        
        % Only set the input current to each unit as the input current when
        % the stimulus is on from hstart to hend.
        if ( t(i) >= hstart ) && (t(i) <= hend )
            IstimE = hE;                        % Set input to E cells
            IstimI = hI;                        % Set input to I cells
        else
            IstimE = zeros(1,Ncells);           % Otherwise no input current
            IstimI = zeros(1,Ncells);           % Otherwise no input current
        end
        IstimE = IstimE + 0.1*randn(1,Ncells);
        IstimI = IstimI + 0.1*randn(1,Ncells);
        % Update rates of all excitatory units based on rates of all units
        % at the previous time point
        rE(i,:) = rE(i-1,:)*(1-dt/tauE) + ...
            dt*(IstimE + rE(i-1,:)*WEE + rI(i-1,:)*WIE + I0E)/tauE;
        
        % Update rates of all inhibitory units based on rates of all units
        % at the previous time point
        rI(i,:) = rI(i-1,:)*(1-dt/tauI) + ...
            dt*(IstimI + rE(i-1,:)*WEI + rI(i-1,:)*WII +I0I)/tauI;
        
        rE(i,:) = max(rE(i,:),0);               % Rates cannot be less than 0
        rI(i,:) = max(rI(i,:),0);               % Rates cannot be less than 0
        
    end
    
    %% Now plot the results for the contrast used
    figure(1)
    if ( trial == 1 )
        clf;                    % Clear figure on first trial
    end
    
    subplot(2,1,1)
    % plot rate of all E cells at end of simulation
    plot(rE(Nt,:),'g')
    hold on
    xlabel('cell index')
    ylabel('rate of E-neurons')
    
    subplot(2,1,2)
    % plot rate of all I cells at end of simulation
    plot(rI(Nt,:),'g')
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
    plot(t,rI(:,cuecell))
    hold on
    % plot rate of I cell with null direction input as a function of time
    plot(t,rI(:,1+mod(cuecell+Ncells/2-1,Ncells)))
    xlabel('time')
    ylabel('rate of inhibitory cell')
    
    figure(3)
    if ( trial == 1 )
        clf;                    % Clear figure on first trial
    end
    
end % Next contrast

figure(23)
switch question_part
    case 'A'
        subplot(1,2,1)
        mesh(WEE)
        colormap(gray)
        axis([0 50 0 50 0 0.32])
        xlabel('Neuron i')
        ylabel('Neuron j')
        zlabel('W^{EE}_{ij}')
    case 'C'
        subplot(1,2,2)
        mesh(WII)
        colormap(gray)
        axis([0 50 0 50 -0.4 0])
        xlabel('Neuron i')
        ylabel('Neuron j')
        zlabel('W^{II}_{ij}')
end
annotation('textbox',[0 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.5 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')


figure(24)
switch question_part
    case 'A'
        subplot('Position',[0.12 0.72 0.85 0.22])
 
    case 'B'
        subplot('Position',[0.12 0.4 0.85 0.22])
        
    case 'C'
        subplot('Position',[0.12 0.08 0.85 0.22])        
end
imagesc(rE')

colormap(gray)
colorbar
title(colorbar,'Rate (Hz)')
ylabel('Unit label')
xlabel('Time (sec)')
set(gca,'XTick',[1 round(hstart/dt) round(hend/dt) 3*Nt/4 Nt])
set(gca,'XTickLabel',{'0' num2str(hstart) num2str(hend) num2str(3*tmax/4) num2str(tmax)})
caxis([0 32])
annotation('textbox',[0 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.0 0.66 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.0 0.34 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')


figure(25)
clf
subplot('Position',[0.08 0.23 0.12 0.6])
plot(WEE(end/2,:),'k')
hold on
plot([0 Ncells], [0 0], 'k:')
axis([ 0 Ncells -0.35 0.35])
set(gca,'XTick',[0 Ncells/2 Ncells])
set(gca,'XTickLabel',{'-\pi/2' '0' '\pi/2'})
xlabel('Tuning Difference')
ylabel('Connection Weight')
title('E-to-E')

subplot('Position',[0.27 0.23 0.12 0.6])
plot(WEI(end/2,:),'k')
hold on
plot([0 Ncells], [0 0], 'k:')
axis([ 0 Ncells -0.35 0.35])
set(gca,'XTick',[0 Ncells/2 Ncells])
set(gca,'XTickLabel',{'-\pi/2' '0' '\pi/2'})
xlabel('Tuning Difference')
title('E-to-I')

subplot('Position',[0.46 0.23 0.12 0.6])
plot(WIE(end/2,:),'k')
hold on
W_di_EE = WEI*WIE;
plot([0 Ncells], [0 0], 'k:')
axis([ 0 Ncells -0.35 0.35])
set(gca,'XTick',[0 Ncells/2 Ncells])
set(gca,'XTickLabel',{'-\pi/2' '0' '\pi/2'})
xlabel('Tuning Difference')
title('I-to-E')

subplot('Position',[0.65 0.23 0.12 0.6])
plot(W_di_EE(end/2,:),'k')
hold on
plot([0 Ncells], [0 0], 'k:')
axis([ 0 Ncells -0.35 0.35])
set(gca,'XTick',[0 Ncells/2 Ncells])
set(gca,'XTickLabel',{'-\pi/2' '0' '\pi/2'})
xlabel('Tuning Difference')
title('Disynaptic E-to-E')

subplot('Position',[0.84 0.23 0.12 0.6])
plot(WEE(end/2,:)+W_di_EE(end/2,:),'k')
hold on
plot([0 Ncells], [0 0], 'k:')
axis([ 0 Ncells -0.35 0.35])
set(gca,'XTick',[0 Ncells/2 Ncells])
set(gca,'XTickLabel',{'-\pi/2' '0' '\pi/2'})
xlabel('Tuning Difference')
title('Total E-to-E')

annotation('textbox',[0 1 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.23 1 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.42 1 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.61 1 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','D')
annotation('textbox',[0.8 1 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','E')
