% Figure_8_7.m
%
% Assume postsynaptic spike with stereotypical membrane potential
% deflection at a time t=0.
% Assess impact on calcium concentration of a presynaptic spike at
% different relative times.
% Impact of the presynaptic spike on the membrane potential is ignored in
% this oversimplified illustration.
%
% This code was used to produce Figure 8.7 in the book
% An Introductory Course in Computational Neuroscience,
% by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
    figure(1)
    clf
%%
% NMDA receptor activation depends on the concentration of magnesium ions,
% since these are the cause of the voltage-dependent block
Mg_conc = 1e-3;
% f_of_v is the voltage-dependence of NMDA opening
f_of_v = @(V) 1./(1+Mg_conc*exp(-V/0.01613)/0.00357 );

tau_Ca = 0.050;             % calcium decay time constant
tau_glutamate = 0.050;      % glutamate unbinding from NMDA time constant

%% Simulation setup.
% All times are relative to the time of the postsynaptic spike at t=0.
pre_times = [-0.005 0.005]; % relative times of presynaptic spikes to be tested
Ntrials = length(pre_times);    % one trial for each pre-post timing

dt = 0.0001;                % time-step for simulation (sec)
tmin = -0.050;              % time to commence simulation
tmax = 0.200;               % time to end simulation
tvec = tmin:dt:tmax;        % vector of time-points
Nt = length(tvec);          % number of time-points

%% Parameters and variables for the neuron
E_L = -0.075;               % leak potential
Vth = -0.055;               % threshold potential
Vreset = -0.085;            % reset voltage
time_factor = 0.5;          % scales time-dependence of V
Vfactor = 0.0001;           % sets scale of membrane potential
% In the next line an approximate trace of membrane potential versus time is
% produced for the postsynaptic cell, with a spike at t=0.
V = (Vreset+Vth)/2 + Vfactor*tan( (log(tvec+1)/time_factor) + pi/2) ;
% Clamp the membrane potentiak between its reset potential and a maximum
V = min(V,0.15);
V = max(V,Vreset);
plot(tvec,V*1000)           % plot the approximated voltage trace

NMDA_open = zeros(Ntrials,Nt);  % fraction of channel opening
Ca = zeros(size(NMDA_open));    % postsynaptic calcium concentration

%% Now loop through trials, simulating the postsynaptic response
for trial = 1:Ntrials
    input_time = pre_times(trial);      % time of presynaptic spike
    % relative amount of glutamate bound to the postsynaptic receptor
    glutamate_bound = (tvec > input_time).*exp( (input_time-tvec)/tau_glutamate);
    % Channel opening is product of the binding of glutamate and a function
    % of postsynaptic membrane potential that determines the amount of
    % magnesium block of the channel
    NMDA_open(trial,:) = glutamate_bound.*f_of_v(V);
    
    % Calcium flows through the channels at a rate proportional to their open
    % probabilities and decays exponentially with a fixed time constant
    for i = 2:length(Ca)
        Ca(trial,i) = Ca(trial,i-1)*exp(-dt/tau_Ca) + NMDA_open(trial,i);
    end
    
    % Now plot all the results
    figure(1)
 
    subplot('Position',[trial*0.5-0.425 0.7 0.35 0.25])
    yyaxis left
    plot(tvec,glutamate_bound,'k:')
    hold on
    axis([-0.015 0.075  0 1])
    ylabel('Glu bound')
    
    ax = gca;
    ax.YColor = 'k';
    yyaxis right
    plot(tvec,V*1000,'k');
    ylim([-85 15])
    ax = gca;
    ax.YColor = 'k';

    legend('Glu bound', 'V_{m}')
    if ( trial == 1 )
        title('Pre-before-post')
    end
    if ( trial == Ntrials )
        title('Post-before-pre')
    end
    
    subplot('Position',[trial*0.5-0.425 0.39 0.35 0.25])
    plot(tvec,NMDA_open(trial,:),'k')
    axis([-0.015 0.075 0 1]);
    hold on
    ylabel('Open fraction')
    
    subplot('Position',[trial*0.5-0.425 0.08 0.35 0.25])
    plot(tvec,Ca(trial,:),'k')
    
    axis([-0.015 0.075  0 15]);
    xlabel('Time (sec)')
    ylabel('[Ca]')
    
end


%% Finally label the figures A1,A2-C1,C2 at the appropriate points
annotation('textbox',[0.00 0.985 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A1')
annotation('textbox',[0.5 0.985 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A2')
annotation('textbox',[0.00 0.66 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B1')
annotation('textbox',[0.5 0.66 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B2')
annotation('textbox',[0.00 0.35 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C1')
annotation('textbox',[0.5 0.35 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C2')
annotation('textbox',[0.425 0.97 0.02 0.03],'LineStyle','none', ...
    'FontSize',14,'FontWeight','Bold','String','V_{m} (mV)')
annotation('textbox',[0.925 0.97 0.02 0.03],'LineStyle','none', ...
    'FontSize',14,'FontWeight','Bold','String','V_{m} (mV)')
