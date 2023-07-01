% chaotic_3units.m
%
% This code is an example of a chaotic circuit produced by three coupled
% firing-rate model units.
% Units have saturating-linear f-I curves.
%
%  A small change in initial conditions causes
%  changes in the behavior that accumulate over time, so that for
%  example the timings of small peaks of activity in r1 eventually
%  become unpredictable. Note the near periodicity of chaotic
%  activity that can produce peaks in the power spectrum.
%
% The code is altered from that found in:
% Dynamical systems, attractors, and neural circuits. P Miller (2016). 
% F1000Res. 2016 May 24;5. pii: F1000 Faculty Rev-992.
%
% This code was used to produce Figure 7.12 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Set up the time vector
tmax = 1.5;             % maximum time to simulate   
dt = 0.0001;            % time-step
tvec = 0:dt:tmax;       % vector of time-points
Nt = length(tvec);      % number of time-points to simulate

%% Network and unit parameters

Nunits = 3;             % Number of units to simulate
r = zeros(Nunits,Nt);   % array of rate of each cell as a function of time
rmax = 100;             % maximum firing rate
tau = 0.010;            % base time constant for changes of rate
% Ithresh is threshold current needed for spiking. 
% Negative values mean spontaneous rate.
Ithresh = [-10; 12; -25 ]; 

% 3x3 connectivity matrix (row = presynaptic, column = postsynaptic unit)
W = [2.05 -1.6 0.0 ; 1.2 0.0 0.2; -0.1 -5.4 2.2 ];

% Two sets of initial conditions with a tiny difference between them
rinit1 = [15.01; 10; 30];
rinit2 = [15; 10; 30];

Ntrials = 2;            % Number of trials (one per initial condition)

%% Figure positions for subplots
figure(1);
subfigB = [0.16 0.4 0.78 0.2];
subfigC = [0.16 0.08 0.78 0.2];
set(0,'DefaultLineLineWidth',2,...      % Plotting parameters
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
rplotmax = 102;                 % sets max of y-axis on r vs t figures

%% Now loop through two trials
for trial = 1:Ntrials
    
    % Set initial conditions for each trial
    if ( trial == 1 )
        r(:,1) = rinit1;
    else
        r(:,1) = rinit2;
    end
        
    % Subplots B and C will both be rate vs time.
    if ( ( trial == 1 ) && ( Ntrials > 1 ) )
        clf
        subplot('position',subfigB)
    else
        subplot('position',subfigC)
    end
    axis([0 tmax 0 rplotmax])
    hold on
    xlabel('Time (sec)')
    ylabel('Firing rate (Hz)')
    
    
    %% Now simulate through time
    for i=2:Nt
        I = W*r(:,i-1);                   % total current to each unit
        newr = r(:,i-1) + dt/tau*(I-Ithresh-r(:,i-1));  % Euler-Mayamara update of rates
        r(:,i) = max(newr,0);                           % rates are not negative
        r(:,i) = min(r(:,i),rmax);                      % rates can not be above rmax
        
        %% Now plot figures rate vs time for each color-coded unit in figure 1,
        %  with a new panel for each trial updating every 1/10 of the
        %  simulation
        
        if ( mod(i,round(Nt/20) ) == 0 )
            if ( Ntrials == 1 )
                subplot('position',subfigB)
                plot(tvec(20:20:i),Iapp(1,20:20:i),'k')
                hold on
                plot(tvec(20:20:i),Iapp(2,20:20:i),'r--')
                plot(tvec(20:20:i),Iapp(3,20:20:i),'g:')
            end
            
            if ( ( trial == 1 ) && ( Ntrials > 1 ) )
                subplot('position',subfigB)
            else
                subplot('position',subfigC)
            end
            plot(tvec(20:20:i),r(1,20:20:i),':','Color',[0.25 0.25 0.25])
            hold on
            plot(tvec(20:20:i),r(2,20:20:i),'-','Color',[0.5 0.5 0.5])
            plot(tvec(20:20:i),r(3,20:20:i),'k','LineWidth',3)
            
            pause(0.05); drawnow
        end
    end
    
    if ( trial == Ntrials )
        legend('r(1)', 'r(2)', 'r(3)')
        
    end
end

% Finally add labels to panels in the figure
annotation('textbox',[0.0 0.95 0.05 0.05],'String','A','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.0 0.62 0.05 0.05],'String','B','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.0 0.3 0.05 0.05],'String','C','LineStyle','none','FontSize',16,'FontWeight','bold')

