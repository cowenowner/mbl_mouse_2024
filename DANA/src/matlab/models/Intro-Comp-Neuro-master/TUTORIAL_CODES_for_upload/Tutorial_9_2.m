% Tutorial_9_2.m
%
% This code produces an inhomogeneous Poisson process, with an initial rate
% and a final rate, with a switch between the two rates at a random point
% in the time interval. Across trials the rates are identical, but the
% change-point varies. The log-likelihood ratio of the change-point being
% at a given point in time versus no change-point is calculated as a
% function of time for each trial. The most likely change-point is then
% compared with the actual change-point.
%
%   This code implements a solution of Tutorial 9.2 of Chapter 9 in the
%   textbook,
%   An Introductory Course in Computational Neuroscience
%   by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

dt = 0.001;                         % Time-step
tmax = 20;                          % Maximum time
t = 0:dt:tmax;                      % Vector of time points
Nt = length(t);                     % Number of time points
r0 = 2;                             % Initial rate (fixed across trials)
rf = 10;                            % Final rate (fixed across trials)
Ntrials = 50;                       % Number of trials

spikes = zeros(Ntrials,Nt);         % Will contain spike train on each trial
r=zeros(Ntrials,Nt);                % Underlying rate on each trial
jumptime = tmax*rand(1,Ntrials);    % Random times of rate-change on each trial
logP = zeros(Ntrials,Nt);           % Log likelihood of change-point at each time per trial
detectjump = zeros(1,Ntrials);      % Time of estimated rate-change on each trial (max of logP)

for trial = 1:Ntrials               % Simulate many trials with different random spike distributions
    display(trial)          % Print to screen to indicate progress                 
    
    for i = 1:length(t);            % Loop through time vector
        if t(i) < jumptime(trial)   % If time is before the change-point
            r(trial,i) = r0;        % rate is initial rate
        else;                       % otherwise
            r(trial,i) = rf;        % rate is final rate
        end
    end
    % Next line produces spikes as a Poisson process with a rate that can
    % vary through the trial
    spikes(trial,:) = rand(1,Nt) < dt*r(trial,:);
        
    % Now calculate just using the spike train, where the rate-change might
    % be.
    for changepoint = 2:Nt-1;   % Loop through all time points
        
        N1 = sum(spikes(trial,1:changepoint));      % No. of spikes before time-point
        N2 = sum(spikes(trial,changepoint+1:end));  % No. of spikes after time-point
        
        T1 = dt*(changepoint);                      % Time before time-point
        T2 = dt*(Nt-changepoint);                   % Time after time-point
            
        r1_opt = N1/T1;                             % Mean rate before time-point
        r2_opt = N2/T2;                             % Mean rate after time-point
        
        % Next line is the log-likelihood of observing these two spike
        % counts at the different rates versus at the constant mean rate
        logP(trial,changepoint)  = N1*log(r1_opt) + N2*log(r2_opt) -N1-N2;     
        
    end
    
    % Find the maximum likelihood
    [val tjump] = max(logP(trial,2:end));
    
    % Record the estimated change-point for that trial
    detectjump(trial) = tjump*dt;

    % Next line displays actual then estimated changepoint
    data = [jumptime(trial) detectjump(trial)]
    
   
    
end

% Next line calculates correlation between actual and estimated
% change-points
correlation = corr(jumptime',detectjump')

% Now plot estimated versus actual change points for all trials
figure()
plot(jumptime,detectjump,'x')
xlabel('Actual rate transition')
ylabel('Estimated rate transition')



