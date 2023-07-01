% single_jumprate_detect.m
%
% A code that generates a single inhomogeneous Poisson spike train, at a
% rate that changes discretely only once between two fixed values. 
%
% Then an algorithm calculates the log likelihood of the jump occurring at
% any point in the simulation.
%
% This code was used to produce Figure 9.5 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;              % remove all prior variables and parameters from memory

dt = 0.001;         % time step for binning spikes and detecting rate-change
tmax = 10;          % total time of spike-train
t = 0:dt:tmax;      % vector of time bins
Nt = length(t);     % number of time bins
r0 = 5;             % initial rate of Poisson process
rf = 10;            % final rate of Poisson process
jumptime = tmax*0.6;    % actual simulated time of rate-change

spikes = zeros(1,Nt);   % will contain a 1 for each spike
r=zeros(1,Nt);          % rate of Poisson process
% logP is log of probability of spike train given rate jump at that time-point
logP = zeros(1,Nt);     
llr = zeros(1,Nt);  % log-likelihood ratio (comparing logP to no rate-jump)

% now simulate the spike train
for i = 1:length(t)
    % determine the Poisson rate at each time-point
    if t(i) < jumptime
        r(i) = r0;
    else
        r(i) = rf;
    end
    
    % produce a spike probabilistically at that rate
    if rand(1) < dt*r(i)        % using a Poisson spike generator
        spikes(i) = 1.0;  % include a spike if the rate is high enough
    end
    
end

%% In the second part of the code use the spike times to try to extract the 
%  time-point when the rate changed (the chang-epoint)
for changepoint = 2:Nt-1; % loop through all possible change-points
    
    N1 = sum(spikes(1:changepoint));        % no. of spikes before change-point
    N2 = sum(spikes(changepoint+1:end));    % no. of spikes after change-point
    N = N1+N2;                              % total no. of spikes
    
    T1 = dt*(changepoint);                  % time before change-point
    T2 = dt*(Nt-changepoint);               % time after change-point
    
    r1_opt = N1/T1;                 % optimal Poisson rate before change-point
    r2_opt = N2/T2;                 % optimal Poisson rate after change-point
    
    % calculate log Prob of spike trains if no. of spikes > 0
    % log prob of spike train before change-point
    if ( r1_opt > 0 )
        logP(changepoint)  = N1*log(r1_opt);
    else
        logP(changepoint)  = 0;
    end
    % add log prob of spike train after change-point
    if ( r2_opt > 0 )
        logP(changepoint)  = logP(changepoint) ...
            + N2*log(r2_opt);
    else
        logP(changepoint)  = logP(changepoint);
    end
    
end
% for log-likelihood ratio, subtract log prob of spike train at fixed rate
llr(2:Nt-1) = logP(2:Nt-1)+N*log(tmax)-N*log(N);

% estimate the time-point of the transition from the maximum of logP
[val tjump] = max(logP(2:end));
detectjump = tjump*dt;

% print to screen actual and estimated change-points
data = [jumptime detectjump] 

%% Finally plot the log-likelihood ratio of a transition at a given 
%  time-point versus no transition as a function of time of change-point
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf
% include the spike times in the final figure 
spiketimes = dt*find(spikes);
subplot('Position',[0.1 0.75 0.85 0.2])
fill([jumptime jumptime tmax tmax], [-0.5 1.5 1.5 -0.5],[0.75, 0.75, 0.75], ...
    'EdgeColor','none')
hold on
set(gca,'Layer','top')

for i = 1:N
    plot([spiketimes(i); spiketimes(i)],[0; 1],'k')
    hold on
end
axis([0 tmax -0.5 1.5])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
ylabel('Spikes')

subplot('Position',[0.1 0.15 0.85 0.6])

% shade the time after the actual transition in gray
fill([jumptime jumptime tmax tmax], [0 ceil(max(llr)) ceil(max(llr)) 0],[0.75, 0.75, 0.75], ...
    'EdgeColor','none')
hold on
set(gca,'Layer','top')

plot(t,llr,'k')
xlabel('Time (sec)')
ylabel('Log likelihood ratio')

axis([0 tmax 0 ceil(max(llr))])

