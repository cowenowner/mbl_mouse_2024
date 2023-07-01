% Figure_6_9.m
%
% Firing rate model for two populations with recurrent excitation,
% coupled by cross-inhibition.
%
% This code was used to generate Figure 6.9 of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

%% Set up the simulation parameters
dt = 0.0005;                % Time step to use
tmax = 3;                   % Maximum time for simulation
tvec = dt:dt:tmax;          % Vector of time points
Nt = length(tvec);          % Number of time points
Ntrials = 100;              % Number of different simulation trials
stim_noise_on_flag = 1;     % If "1" then noise is in the stimulus
V_noise_on_flag = 0;        % If "1" then additional noise in each unit

%% Parameters for firing-rate response curve which is the model by
%  FS Chance and LF Abbott, Prog Brain Res 149:147-155, 2005
%  The model approximates well the output of a leaky integrate-and-fire
%  (LIF) neuron given its steady state voltage. So LIF parameters are
%  needed.
delta_V = 0.025;            % Threshold voltage minus reset voltage
% The parameters a and sigma_V both set the curvature of firing rate curve.
% sigma_V should be the standard deviation of voltage noise and then gets
% scaled (inversely) by a. Here sigma_V does not actually add noise
% fluctuations unless the flag V_noise_on_flag is set to 1.
a = 0.05;
sigma_V = 0.0002;
Vth = -0.050;               % Threshold voltage (V)
E_L = -0.065;               % Leak voltage (V)
g_L = 0.06;                 % Leak conductance (nS)
E_I = -0.065;               % Inhibitory synapse reversal potential (V)
tau = 0.010;                % Time constant for firing-rate to change
decision_threshold = 40;    % Threshold indicates decision is made

if ( V_noise_on_flag )
    V_Noise = sigma_V;
else
    V_Noise = 0;                % Set voltage noise to zero
end

g_stim = 0.03;

%% Parameters that vary across conditions are given in vectors, with each
%  element of the vector to be used on one corresponding trial.
%  First trial is depression only.
%  Second trial is facilitation only.
%  Third trial has both facilitation and depression.


%% Synaptic response parameters

alpha0 = 0.5;      % Baseline postsynaptic amplitude.
Wee = 0.8;          % Strength of synaptic feedback.
Wei = 8.5;        % Strength of synaptic feedback.
f_f = 0;          % Facilitation factor (0 means no facilitation)
p0 = 0.05;        % Baseline vesicle release probability.

% Time constant for synaptic depression, zero means no depression
tau_d = 0.1;

tau_sE = 0.050;                 % Synaptic time constant for E-feedback
tau_sI = 0.005;                 % Synaptic time constant for I-feedback
tau_f = 0.3;                    % Facilitation time constant

%% Stimulus properties that can vary per trial
stim_amp = 0.25;                % Mean stimulus amplitude
stim_diff = 0.02;               % Amplitude difference between 2 stimuli
t_on = 0.5;                     % Time of stimulus onset
stim_length = tmax-t_on;        % Duration of stimulus
t_off = t_on + stim_length;     % Time for stimulus offset

if ( stim_noise_on_flag )           % If the stimulus is noisy
    stim_noise_amp = 0.01;         % Amplitude of stimulus noise
else
    stim_noise_amp = 0;             % Otherwise set to zero
end

% Because the stimulus noise scales to large values with small dt, it is
% computationally safer to keep the noise timescale fixed at a larger value
% than dt.
noise_time_scale = 0.005;           % 5ms is smaller than most cellular timescales
% noise_bin is the number of time steps in the noise timescale, for which
% the noise input will be identical
noise_bin = round(noise_time_scale/dt);

% N_noise_bins is number of distinct noise values
N_noise_bins = ceil(stim_length/noise_time_scale);

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
% LineStyle_vec is to allow different linestyles within a loop
% Note { } is used for a vector of character strings
LineStyle_vec = {'-', '--', ':'};
rmax = decision_threshold + 10;     % Maximum rate for the plots

Wee_mat = [Wee 0; 0 Wee];   % Matrix of E-to-E synaptic strengths
Wei_mat = [0 Wei; Wei 0];   % Matrix of E-to-I synaptic strengths

% Set the stimuli such that stim_amp is the mean and the difference
% (stim_amp1 - stim_amp2) is obtained from stim_diff_vec
% stim_amp1 is amplitude of stimulus to Unit 1
stim_amp1 = stim_amp + 0.5*stim_diff;
% stim_amp2 is amplitude of stimulus to Unit 2
stim_amp2 = stim_amp - 0.5*stim_diff;

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

% The variable "counts" is to record number of times each Unit is the
% "winner".
counts = zeros(1,2);        

% Loop through trials of identical stimuli to accumulate statistics
for trial = 1:Ntrials;
    stim = zeros(2, Nt);        % Initialize vector for the stimulus
    % then add stimulus to appropriate time points.
    stim(1,round(t_on/dt):round(t_off/dt)) = stim_amp1;
    stim(2,round(t_on/dt):round(t_off/dt)) = stim_amp2;
    
    % Now add noise to the stimulus while it is present
    % Loop through all the separate noise bins
    for bin = 1:N_noise_bins;
        % i_on and i_off indicate the first and last time-points in
        % a bin with fixed noise
        i_on = round(t_on/dt) + (bin-1)*noise_bin+1;
        i_off = round(t_on/dt) + bin*noise_bin;
        % Ensure a separate noise-value is given to each unit but
        % that value is fixed across the set of time-bins within
        % one noise-bin. The value is scaled by the amplitude of
        % noise and inversely by the square-root of the noise-bin
        % size.
        stim(:,i_on:i_off) = stim(:,i_on:i_off) + ...
            stim_noise_amp*randn(2,1)*ones(1,noise_bin)/sqrt(noise_time_scale);
    end
    
    % Below is the formula for the firing rate response function that will
    % be used in both the steady state curves and the simulation.
    % This is a function of "x" where "x" is the value of the
    % variable (or vector of values) sent to the function.
    r_of_Vss = @(x) (x - Vth)./(tau*delta_V*(1-exp(-a*(x-Vth)/sigma_V)));
    
    r = zeros(2,Nt);        % Initialize vector for firing rates
    sE = zeros(2,Nt);       % Initialize vector for synaptic inputs
    sI = zeros(2,Nt);       % Initialize vector for synaptic inputs
    D = ones(2,Nt);         % Initialize vector for depression variable
    F = ones(2,Nt);         % Initialize vector for facilitation variable
    
    % Now produce a vector of stimulus-independent noise for each cell
    V_Noise_vec = V_Noise*randn(2,Nt)*sqrt(tau/dt);
    
    %% Now carry simulate the same model as a function of time
    % In models with multiple variables, each variable is updated using
    % the values of the set of variables on the previous time step
    
    i = 1;              % i denotes the time-step
    % A while loop is used to exit the trial once threshold is
    % reached
    while ( ( i < Nt) && (max(r(:,i)) < decision_threshold ) )                               % Loop through time points
        i = i+1;
        
        % Calculate the steady-state membrane potential for a
        % leaky-integrate-and-fire cell based on the excitatory and
        % inhibitory synaptic inputs.
        % Vss is a column vector with one value for each cell.
        % The synaptic inputs are obtained by matrix multiplication
        % using the connectivity matrices Wee_mat and Wei_mat, plus
        % the stimulus input conductance in "stim".
        Vss = (E_L*g_L +E_I*Wei_mat*sI(:,i-1) )./ ...
            (g_L  + Wee_mat*sE(:,i-1) + Wei_mat*sI(:,i-1) +g_stim*stim(:,i));
        
        % Now update the firing rate using the equation
        % dr/dt = -(r + r(Vss) )/tau
        r(:,i) = r(:,i-1)*(1-dt/tau) + ...
            dt/tau*r_of_Vss(Vss+V_Noise_vec(:,i));
        
        % Next update the facilitation variable which depends on firing
        % rate alone, using the equation
        % dF/dt = (1-F)/tau_F + f_f*(1/p0 - F)*r
        % Note that if f_f is zero then F is fixed at 1.
        F(:,i) = F(:,i-1) +  dt/tau_f*(1-F(:,i-1)) + ...
            f_f*(1/p0-F(:,i-1))*dt.*r(:,i-1);
        
        % Next update the depression variable, which depends on both firing
        % rate and the facilitation variable, using the equation
        % dD/dt = (1-D)/tau_D - pr.D.r (where pr = p0.F)
        if ( tau_d > 0 )        % tau_d = 0 ensures D = 1 always
            D(:,i) = D(:,i-1) +(1-D(:,i-1))*dt/tau_d - ...
                p0*F(:,i-1).*D(:,i-1).*r(:,i-1)*dt;
        else
            D(:,i) = 1;     % No impact of depression if D = 1.
        end
        
        % Finally update the synaptic response variable,s, which depends on
        % all of rate, facilitation and depression variables, using the
        % equation
        % ds/dt = -s/tau_s + alpha.pr.r.(1-s)
        % where pr = p0.F and alpha = alpha_0*D
        sE(:,i) = sE(:,i-1)*(1-dt/tau_sE) + ...
            dt*F(:,i-1)*p0.*D(:,i-1)*alpha0.*r(:,i-1) ...
            .*(1-sE(:,i-1));
        % The inhibitory synapti response differs because of the
        % different time constant. Here, the presynaptic properties
        % (F and D) are identical for inhibitory and excitatory
        % synapses.
        sI(:,i) = sI(:,i-1)*(1-dt/tau_sI) + ...
            dt*F(:,i-1)*p0.*D(:,i-1)*alpha0.*r(:,i-1) ...
            .*(1-sI(:,i-1));
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
        decision_time(trial) = dt*time_point-t_on;
    end
    
    %% Now plot the inidividual trial results
    
    if ( ( trial < 5 ) )
        if ( trial == 1 )
            figure(1)       
            clf
            subplot('Position',[0.6 0.6 0.37 0.37])

            % produce a solid gray rectangle across the region of the stimulus
            fill([0 0 t_off-t_on t_off-t_on],[0 rmax rmax 0],[0.75, 0.75, 0.75], ...
                'EdgeColor','none')
            hold on;        % do not overwrite the gray rectangle
            
            subplot('Position',[0.1 0.6 0.37 0.37])
           % produce a solid gray rectangle across the region of the stimulus
            fill([0 0 t_off-t_on t_off-t_on],[0 rmax rmax 0],[0.75, 0.75, 0.75], ...
                'EdgeColor','none')
            hold on
            plot(tvec-t_on,r(1,:),'k');     % plot firing rate as a function of time
            plot(tvec-t_on,r(2,:),'k:');    % plot firing rate as a function of time
            hold on
            xlabel('Time (sec)')            % Label x-axis on last trial
            ylabel('Firing Rate (Hz)')
            axis([-0.25 1.25 0 45])         % Set the range of each axis
            set(gca,'Layer','top')          % Ensures axes lines are visible on top of the plot



        end
        figure(1)
        subplot('Position',[0.6 0.6 0.37 0.37])
        plot(tvec-t_on,r(1,:),'k');         % plot firing rate as a function of time
        xlabel('Time (sec)')                % Label x-axis on last trial
        ylabel('Firing Rate (Hz)')
        axis([-0.25 1.25 0 45])             % Set the range of each axis
        set(gca,'Layer','top')              % Ensures axes lines are visible on top of the plot

        
        figure(2)
        plot(tvec-t_on,r(2,:),'k');         % plot firing rate as a function of time
        hold on
        xlabel('Time (sec)')                % Label x-axis on last trial
        ylabel('Firing Rate (Hz)')
        axis([-0.25 1.25 0 rmax])           % Set the range of each axis
        set(gca,'Layer','top')              % Ensures axes lines are visible on top of the plot

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
end;            % Go to next trial

% max_i is the greatest time-point to reach threshold across all trials
% that Unit 1 was the "winner".
max_i = max(find(trial_counts));
% Now divide the sum of rates as a function of time, across trials by the 
% number of trials still active at that point in time. The last index is
% needed to prevent division by zero beyond the slowest decision time.
mean_r(1,1:max_i) = mean_r(1,1:max_i)./trial_counts(1:max_i);
mean_r(2,1:max_i) = mean_r(2,1:max_i)./trial_counts(1:max_i);

% For the reverse rates, all trials are active -- the high indices
% correspond to spontaneous activity before stimulus, which we do not care
% about, so is set to zero.
mean_reverse_r = mean_reverse_r./counts(1);

% First plot stimulus-aligned firing rates of the two units
figure(1)
subplot('Position',[0.1 0.1 0.37 0.37])
% The next line gives the shaded gray box a handle called "shade" that can
% then be used to switch off a legend entry for the object
shade = fill([0 0 t_off-t_on t_off-t_on],[0 rmax rmax 0],[0.75, 0.75, 0.75], ...
    'EdgeColor','none')
set(get(get(shade,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude shaded box from legend
hold on
plot(tvec-t_on,mean_r(1,:),'k');
plot(tvec-t_on,mean_r(2,:),'k:');
xlabel('Time from stimulus (sec)')
ylabel('Mean Rate (Hz)')
axis([-0.25 1.25 0 45])
set(gca,'Layer','top')  % Ensures axes lines are visible on top of the plot
legend('Unit 1','Unit 2')

% Next plot response-aligned firing rates by using the reverse-time rates
% in order of descending index. Index 1 corresponds to the
% threshold-crossing of each trial and higher indices are preceding
% time-points.
figure(1)
subplot('Position',[0.6 0.1 0.37 0.37])
fill([ -t_off+t_on -t_off+t_on 0 0],[0 rmax rmax 0],[0.75, 0.75, 0.75], ...
    'EdgeColor','none')
hold on
plot(tvec-tmax,mean_reverse_r(1,end:-1:1),'k');
plot(tvec-tmax,mean_reverse_r(2,end:-1:1),'k:');
set(gca,'Layer','bottom')  % Ensures axes lines are visible on top of the plot

xlabel('Time before threshold (sec)')
ylabel('Mean Rate (Hz)')
axis([-1.5 0 0 45])
set(gca,'Layer','top')  % Ensures axes lines are visible on top of the plot

annotation('textbox',[0 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.5 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0 0.48 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.5 0.48 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','D')
