% Figure_3_7.m
%
% This code produces five different spike trains with different statistics
% then generates the distribution of inter-spike intervals (ISIs) to be
% plotted and calculates the coefficient of variation (CV) of the ISI 
% distribution and CV_2 (mean coefficient of variation between neighboring
% ISIs).
%
% This code was used to produce Figure 3.7 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initially set the default plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf
%% Now set up common details of the spike trains
clear

Ntrials = 5;                % Number of spike trains to produce
dt = 0.001;                 % Width of time-bin is 1ms
tmax = 300;                 % 300-secs of spike train to generate sufficient statistics
tvec = dt:dt:tmax;          % Vector of time points
Nt = length(tvec);          % NUmber of time points

spikes = zeros(Ntrials,Nt); % Vector of spike times, one row per spike train

% Trial 1 has regular spikes every 50ms
spikes(1,25:50:Nt) = 1;

% Trial 2 has spikes with some variation, sigma about a mean of 50ms
sigma = 0.02;
spike_time = dt;
while spike_time < tmax
    spikes(2,floor(spike_time/dt)) = 1;
    spike_time = spike_time + max(0.05 + sigma*randn(),dt);
end

% Trial 3 is a Poisson distribution (independently random every time-bin)
spikes(3,:) = rand(1,Nt) < dt*20;

% Trial 4 has a step change in spike rate (no noise)
spike_time = 0;
while spike_time < 3*tmax/4
    spike_time = spike_time + 0.15;
    spikes(4,floor(spike_time/dt)) = 1;
end
while spike_time < tmax-0.01
    spike_time = spike_time + 0.01;
    spikes(4,floor(spike_time/dt)) = 1;
end

% Trial 5 is a bursting neuron, with regularly spaced triplets
spikes(5,21:150:Nt) = 1;
spikes(5,25:150:Nt) = 1;
spikes(5,30:150:Nt) = 1;

%% Now analyze each spike train in turn
for trial = 1:Ntrials
    
    spiketimes=dt*find(spikes(trial,:));    % A list of spike times
    isis=1000*diff(spiketimes);             % ISIs converted to ms
    cv=std(isis,1)/mean(isis)               % CV is standard deviation divided by mean
    
    % Now calculate the CV of each consecutive pair of ISIs and sum them
    sum2=0.0;
    for i=2:length(isis)
        sum2 = sum2 + sqrt((isis(i)-isis(i-1))*(isis(i)-isis(i-1))) ...
            /(isis(i)+isis(i-1))*(2.0);
    end
    cv2 = sum2 /(length(isis)-1) % Divide by number of such pairs to get the mean  
    
    % Plot the spike train in a small window
    subplot('Position',[0.05 1.05-trial*0.195  0.55 0.09 ])
    plot(1000*tvec,spikes(trial,:),'k')
    if ( trial == 4 )
        axis ([224250 225250 -0.5 1.5])
        set(gca,'XTick',[224250:200:225250])
    else
        axis([0 1000 -0.5 1.5])
    end
    if ( trial == Ntrials )
        xlabel('Time (ms)')
    end
    set(gca,'XTickLabel',{0 200 400 600 800 1000})
    set(gca,'YTick',[])
    
    switch trial
        case 1
            title('regular')
        case 2
            title('regular + jitter')
        case 3
            title('Poisson')
        case 4
            title('rate jump')
        case 5
            title('3-spike bursts')
    end
    
    % Plot the distribution of ISIs
    subplot('Position',[0.65 1.05-trial*0.195 0.3 0.09 ])
    h=histogram(isis,[5:10:155])
    h.FaceColor = 'k'
    
    axis([0 170 0 2000])
    
    if ( trial == Ntrials )
        xlabel('ISI (ms)')
    end
    % Produce one decimal place (could use format when printing).
    cv = round(cv*10)/10;
    cv2 = round(cv2*10)/10;
    set(gca,'YTick',[])
    title(strcat('   CV=',num2str(cv),', CV_{2}=',num2str(cv2)));
end

%% Label the subplots A-E
annotation('textbox',[0.00 0.985 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A1')
annotation('textbox',[0.00 0.79 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B1')
annotation('textbox',[0.00 0.595 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C1')
annotation('textbox',[0.00 0.40 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','D1')
annotation('textbox',[0.00 0.205 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','E1')
annotation('textbox',[0.60 0.985 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A2')
annotation('textbox',[0.60 0.79 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B2')
annotation('textbox',[0.60 0.595 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C2')
annotation('textbox',[0.60 0.40 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','D2')
annotation('textbox',[0.60 0.205 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','E2')
