% Tutorial_6_2.m
%
% Firing rate model for two populations with recurrent excitation,
% coupled by cross-inhibition to produce a winner-takes-all decision-making 
% network. Uses linear firing-rate curves.
%
% Two modes of operation are possible.
% 1. An integrator mode, in which the rate difference between the two
% groups is proportional to the time-integral of the stimulus difference,
% within bounds (of 0 and rmax).
% 2. A jumping mode, which relies on noise for a transition from a
% low-firing rate state for both units, to a state with one of the two
% units active (a decision state).
%
% This code is a solution of Tutorial 6.2 of the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University (2017)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
mode_number = 2;            % Use "1" for integrator mode or "2" for jumping mode
all_figures_on = 1;

%% Set up the simulation parameters
dt = 0.0005;                    % Time step to use
tmax = 20;                     % Maximum time for simulation
tvec = 0:dt:tmax;              % Vector of time points
Nt = length(tvec);              % Number of time points
Ntrials = 200;                    % Number of different simulation trials
stim_noise_on_flag = 1;         % If "1" then noise is in the stimulus
internal_noise_on_flag = 1;     % If "1" then additional noise in each unit

tau = 0.010;                    % Time constant for firing-rate to change
decision_threshold = 50;        % Threshold indicates decision is made

if ( internal_noise_on_flag )   % evaluates as true if internal_noise_flag is non-zero
    internal_noise = 0.25;
else
    internal_noise = 0;                % Set voltage noise to zero
end

%g_stim_vec = [0.5 1];
g_stim_vec = 0.5
%% Parameters that vary across conditions are given in vectors, with each
%  element of the vector to be used on one corresponding trial.
%  First trial is depression only.
%  Second trial is facilitation only.
%  Third trial has both facilitation and depression.


%% Synaptic response parameters
switch mode_number
    case 1
        Ws = 0.975;           % Strength of self-excitation.
        Wx = -0.025;          % Strength of cross-inhibition.
        Ith = -0.5;
% Results when Ds = 0.5 and s = 1:
% mean_times = 0.5429
% fraction_choose_1 = 0.8152
    case 2
        Ws = 1.05;          % Strength of self-excitation.
        Wx = -0.05;          % Strength of cross-inhibition.
        Ith = 4;
        g_stim_vec = [2 2.5];
% Results when Ds = 0.5 and s = 1:
% mean_times =  0.5409
% fraction_choose_1 = 0.8176
end
N_strengths = length(g_stim_vec);

rmax = decision_threshold + 10;     % Maximum rate for the plots
r_init = -Ith/(1-(Ws+Wx));
r_init = max(r_init,0);
if ( Ws+Wx > 1 ) 
    r_init = 0;
end

%% Stimulus properties that vary per trial

t_on = 1;                           % Time of stimulus onset
stim_amp = 1;                       % Mean stimulus amplitude
stim_diff_vec = [0:0.25:1];      % Amplitude difference between 2 stimuli
stim_length = tmax-t_on;            % Duration of stimulus
t_off = t_on + stim_length;         % Time for stimulus offset

if ( stim_noise_on_flag )           % If the stimulus is noisy
    stim_noise_amp = 0.25;             % Amplitude of stimulus noise
else
    stim_noise_amp = 0;             % Otherwise set to zero
end

% Because the stimulus noise scales to large values with small dt, it is
% computationally safer to keep the noise timescale fixed at a larger value
% than dt.
noise_time_scale = 0.002;           % 2ms is smaller than most cellular timescales
% noise_bin is the number of time steps in the noise timescale, for which
% the noise input will be identical
noise_bin = round(noise_time_scale/dt);

% N_noise_bins is number of distinct noise values
N_stim_noise_bins = ceil(stim_length/noise_time_scale); % Noise in stimulus
N_int_noise_bins = ceil(tmax/noise_time_scale);         % Internal noise

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
% LineStyle_vec is to allow different linestyles within a loop
% Note { } is used for a vector of character strings
LineStyle_vec = {'-', '--', ':'};

%% Set up variables to record results
% Nparameters is the number of different parameter sets where each
% parameter set here contains one value for the difference in stimulus
% amplitudes that corresponds to coherence.
Nparameters = length(stim_diff_vec);

% mean_times will store the mean times to threshold
mean_times = zeros(N_strengths,Nparameters);
% fraction_correct is actually the fraction of trial where Unit 1 reaches
% threshold first
fraction_choose_1 = zeros(N_strengths,Nparameters);
% fraction_undecided is the fraction of trials in which neither Unit
% reaches threshold during the simulation. Ideally this should be zero.
fraction_undecided = zeros(N_strengths,Nparameters);


%% Now loop through different values of the stimulus input conductance that
%  scales the signal and the noise in the stimulus together. A large value
%  of the conductance produces faster responses with more noise, while a
%  low value produces slower responses with less noise impact.
for strength_index = 1:N_strengths
    g_stim = g_stim_vec(strength_index);    % Input conductance
    for parameter = 1:Nparameters           % Loop through stimulus values
        parameter          % Prints to screen to indicate progress
        
        sum_times = 0;          % Initialize for this stimulus set
        % The variable "counts" accumulates the number of times each unit crosses
        % threshold first
        counts = zeros(1,2);    % Initialize to zero for this stimulus set
        
        W_mat = [Ws Wx; Wx Ws];   % Matrix of E-to-E synaptic strengths
        
        % Set the stimuli such that stim_amp is the mean and the difference
        % (stim_amp1 - stim_amp2) is obtained from stim_diff_vec
        % stim_amp1 is amplitude of stimulus to Unit 1
        stim_amp1 = stim_amp + 0.5*stim_diff_vec(parameter);
        % stim_amp2 is amplitude of stimulus to Unit 2
        stim_amp2 = stim_amp - 0.5*stim_diff_vec(parameter);
        
        mean_r = zeros(2,Nt);   % To accumulate trial-averaged rate
        % To accumulate times to threshold, trial-averaged rates will be
        % stored in reverse time order.
        mean_reverse_r = zeros(2,Nt);
        % To count trials that reach each time-point for later
        % normalization to find the mean rates since simulations will stop
        % when a decision is made and the total number of trials is no
        % longer the correct normalization.
        trial_counts = zeros(1,Nt);
        decision_time = zeros(1,Ntrials);   % To store times to threshold.
        
        % Loop through trials of identical stimuli to accumulate statistics
        for trial = 1:Ntrials;
            if ( mod(trial,100) == 0 )
                trial
            end
            
            stim = zeros(2, Nt);        % Initialize vector for the stimulus
            internal_noise_vec = zeros(2, Nt); 
            % then add stimulus to appropriate time points.
            stim(1,round(t_on/dt)+1:round(t_off/dt)) = stim_amp1;
            stim(2,round(t_on/dt)+1:round(t_off/dt)) = stim_amp2;
            
            % Now add noise to the stimulus while it is present
            % Loop through all the separate noise bins
            for bin = 1:N_stim_noise_bins;
                % i_on and i_off indicate the first and last time-points in
                % a bin with fixed noise
                i_on = round(t_on/dt) + (bin-1)*noise_bin+1;
                i_off = min(round(t_on/dt) + bin*noise_bin,Nt);
                % Ensure a separate noise-value is given to each unit but
                % that value is fixed across the set of time-bins within
                % one noise-bin. The value is scaled by the amplitude of
                % noise and inversely by the square-root of the noise-bin
                % size.
                stim(:,i_on:i_off) = stim(:,i_on:i_off) + ...
                    stim_noise_amp*randn(2,1)*ones(1,noise_bin)/sqrt(noise_time_scale);
            end
            
            % Then add internal noise for all time bins
            % Loop through all the separate noise bins
            for bin = 1:N_int_noise_bins;
                % i_on and i_off indicate the first and last time-points in
                % a bin with fixed noise
                i_on = (bin-1)*noise_bin+1;
                i_off = min(bin*noise_bin,Nt);
                % Ensure a separate noise-value is given to each unit but
                % that value is fixed across the set of time-bins within
                % one noise-bin. The value is scaled by the amplitude of
                % noise and inversely by the square-root of the noise-bin
                % size.
                internal_noise_vec(:,i_on:i_off) = internal_noise* ...
                    randn(2,1)*ones(1,noise_bin)/sqrt(noise_time_scale);
            end
            
            
            r = zeros(2,Nt);        % Define vector for firing rates
            r(:,1) = r_init;        % Initialize vector for firing rates
            I = zeros(2,Nt);        % Initialize vector for current inputs
                        
            %% Now simulate the same model as a function of time
            % In models with multiple variables, each variable is updated using
            % the values of the set of variables on the previous time step
            
            i = 1;              % i denotes the time-step
            % A while loop is used to exit the trial once threshold is
            % reached
            while ( ( i < Nt) && (max(r(:,i)) < decision_threshold ) )                               % Loop through time points
                i = i+1;
                I =  W_mat*r(:,i-1) + g_stim*stim(:,i) + internal_noise_vec(:,i);
                
                % Now update the firing rate using the equation
                % dr/dt = -(r + I - Ith )/tau
                r(:,i) = r(:,i-1)*(1-dt/tau) + ...
                    dt/tau*(I-Ith);
                
                r(:,i) = min(r(:,i),rmax);
                r(:,i) = max(r(:,i),0);
                
            end;        % go to next time point
            
            % row and column provided for all values of firing rate above
            % threshold in the next line.
            % the first index "cell_series" is the row label, so is the
            % Unit above threshold
            % the second index "time_series" is the column label so
            % contains the corresponding time points above threshold.
            [cell_series, time_series ] = find(r > decision_threshold);
            
            % in the next line "time_point" is the lowest value of all the
            % time points above threshold and "index" corresponds to which
            % element in the vector "time_series" that it corresponds to.
            % Normally index will be 1.
            [time_point, index] = min(time_series);
            % The Unit that crossed threshold first is "winner"
            winner = cell_series(index);
            
            if ( winner )           % Proceed if there is a decision
                % Update the tally for the winning Unit
                counts(winner) = counts(winner) + 1;
                
                % Update tally of decision-making times
                sum_times = sum_times + dt*time_point-t_on;
                
            end
            
            %% Now plot the results
            if ( all_figures_on )
                if ( ( trial < 5 ))
                    if ( trial == 1 )
                        figure(parameter)       % New figure for each stimulus
                        
                        % produce a solid gray rectangle across the region of the stimulus
                        fill([t_on t_on t_off t_off],[0 rmax rmax 0],[0.75, 0.75, 0.75], ...
                            'EdgeColor','none')
                        hold on;        % do not overwrite the gray rectangle
                        figure(parameter+20)
                        
                    end
                    figure(parameter)
                    plot(tvec,r(1,:),'k');  % plot firing rate as a function of time
                    hold on
                    axis([0 5 0 rmax])       % Set the range of each axis

                    figure(parameter+20)
                    plot(tvec,r(2,:),'k');  % plot firing rate as a function of time
                    hold on
                    xlabel('Time (sec)')        % Label x-axis on last trial
                    axis([0 5 0 rmax])       % Set the range of each axis
                    drawnow;
                end
            end
            % To calculate mean rate we first accumulate across trials then
            % we will noramlize by the number of trials still running at
            % any time point. To achieve the latter we accumulate by 1s in
            % each time-point until the last simulated ("i") in that trial.
            if ( winner == 1 )
                mean_r = mean_r + r;                % accumulate rates
                trial_counts(1:i) = trial_counts(1:i) + 1;  % "i" is last time-point
                % Now store rates in reverse order with the threshold-crossing
                % rate as the initial value, in order to later align to the
                % threshold point.
                mean_reverse_r(:,1:i) = mean_reverse_r(:,1:i) + r(:,i:-1:1);
            end
        end;            % Go to next parameter set and produce a new trial
        mean_times(strength_index, parameter) = sum_times/sum(counts);
        fraction_choose_1(strength_index, parameter)  = counts(1)/Ntrials;
        fraction_undecided(strength_index, parameter)  = 1-sum(counts)/Ntrials;        

        max_i = max(find(trial_counts));
        mean_r(1,1:max_i) = mean_r(1,1:max_i)./trial_counts(1:max_i);
        mean_r(2,1:max_i) = mean_r(2,1:max_i)./trial_counts(1:max_i);
        mean_reverse_r = mean_reverse_r./counts(1);
        
        figure(200)
        hold on
        plot(tvec,mean_r(1,:),'k');
        plot(tvec,mean_r(2,:),'k:');
        figure(201)
        hold on
        plot(tvec-tmax,mean_reverse_r(1,end:-1:1),'k');
        plot(tvec-tmax,mean_reverse_r(2,end:-1:1),'k:');

    end
    
end

if(Nparameters > 1 )
    figure(99)
    subplot(2,1,1)
    legend_vec = [];
    for i = 1:N_strengths
        plot(stim_diff_vec,fraction_choose_1(i,:) ...
            +0.5*fraction_undecided(i,:),'r','linestyle',LineStyle_vec{i});
        hold on
        legend_string = strcat('g_s = ',num2str(g_stim_vec(i)));
        legend_vec{i} = legend_string;
    end
    
    xlabel('Stimulus difference')
    ylabel('Fraction Choose 1')
    legend(legend_vec);
    subplot(2,1,2)
    for i = 1:N_strengths
        plot(stim_diff_vec,mean_times(i,:),'r','linestyle',LineStyle_vec{i});
        hold on
    end
    xlabel('Stimulus difference')
    ylabel('Mean Response Time')
    legend(legend_vec);
    
end
