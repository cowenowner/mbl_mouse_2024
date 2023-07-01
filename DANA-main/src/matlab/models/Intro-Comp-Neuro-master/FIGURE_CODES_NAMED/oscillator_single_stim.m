% oscillator_single_stim.m
%
% Firing rate model for two populations, one excitatory and one inhibitory,
% coupled to produce a PING oscillator with a firing-rate model.
%
% The code calculates the oscillation frequency using two methods, one via
% a Fourier transform. 
%
% This code was used to generate Figure 6.12 of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%% Set up the simulation parameters
dt = 0.0001;                 % Time step to use
tmax = 1.5;                   % Maximum time for simulation
tvec = dt:dt:tmax;           % Vector of time points
Nt = length(tvec);           % Number of time points
N_Units = 2;
Ntrials = 1;               % Number of different simulation trials
stim_noise_on_flag = 0;      % If "1" then noise is in the stimulus
V_noise_on_flag = 0;         % If "1" then additional noise in each unit

%% Parameters for firing-rate response curve which is the model by
%  FS Chance and LF Abbott, Prog Brain Res 149:147-155, 2005
%  The model approximates well the output of a leaky integrate-and-fire
%  (LIF) neuron given its steady state voltage. So LIF parameters are
%  needed.
delta_V = 0.030;            % Threshold voltage minus reset voltage
% The parameters a and sigma_V both set the curvature of firing rate curve.
% sigma_V should be the standard deviation of voltage noise and then gets
% scaled (inversely) by a. Here sigma_V does not actually add noise
% fluctuations unless the flag V_noise_on_flag is set to 1.
a = 1;
sigma_V = 0.001;
Vth = -0.050;               % Threshold voltage (V)
E_L = -0.070;               % Leak voltage (V)
g_L = 0.05;                 % Leak conductance (nS)
E_I = -0.065;               % Inhibitory synapse reversal potential (V)
tau = 0.003;                % Time constant for firing-rate to change

if ( V_noise_on_flag )
    V_Noise = sigma_V;
else
    V_Noise = 0;                % Set voltage noise to zero
end

g_stim = 1;

W_vary = 0;
%% Parameters that vary across conditions are given in vectors, with each
%  element of the vector to be used on one corresponding trial.
%  First trial is depression only.
%  Second trial is facilitation only.
%  Third trial has both facilitation and depression.


%% Synaptic response parameters

alpha0 = 0.5;           % Baseline postsynaptic amplitude.
Wee0 = 25;              % Strength of E-to-E synaptic feedback.
Wei0 = 25;              % Strength of E-to-I synaptic feedback.
Wie0 = 400;             % Strength of I-to-E synaptic feedback.
Wii = 0;                % Strength of I-to-I synaptic feedback.
f_f = 0;                % Facilitation factor (0 means no facilitation)
p0 = 0.1;               % Baseline vesicle release probability.

W_ex_vec = [0.5:0.5:20];

W_in_vec = ones(size(W_ex_vec));

% Time constant for synaptic depression, zero means no depression
tau_d = 0.1;

tau_sE = 0.002;                 % Synaptic time constant for E-feedback
tau_sI = 0.005;                 % Synaptic time constant for I-feedback
tau_f = 0.3;                    % Facilitation time constant

%% Stimulus properties that vary per trial

t_on = dt;                          % Time of stimulus onset
if ( W_vary )
    stim_amp_vec = 0.1*ones(size(W_ex_vec));
else
    stim_amp_vec = [0.1];   % Mean stimulus amplitude
end
stim_length = tmax-t_on;            % Duration of stimulus

if ( stim_noise_on_flag )           % If the stimulus is noisy
    stim_noise_amp = 0.005;         % Amplitude of stimulus noise
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
plot_all = 0;               % Set to 1 to plot each simulation
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
if ( W_vary )
    Nparameters = length(W_ex_vec);
else
    Nparameters = length(stim_amp_vec);
end

mean_f = zeros(1,Nparameters);
Num_cycles = zeros(1,Nparameters);
r_amp1 = zeros(1,Nparameters);
r_amp2 = zeros(1,Nparameters);
mean_r1 = zeros(1,Nparameters);
mean_r2 = zeros(1,Nparameters);

for parameter = 1:Nparameters           % Loop through stimulus values
    parameter          % Prints to screen to indicate progress
    
    sum_times = 0;          % Initialize for this stimulus set
    sum_sqr_times = 0;      % Initialize for this stimulus set
    % The variable "counts" accumulates the number of times each unit crosses
    % threshold first
    counts = zeros(1,2);    % Initialize to zero for this stimulus set
    if ( W_vary)
        Wee = Wee0*W_ex_vec(parameter);
        Wei = Wei0*W_ex_vec(parameter);
        Wie = Wie0*W_in_vec(parameter);
    else
        Wee = Wee0;
        Wei = Wei0;
        Wie = Wie0;
    end
    Wee_mat = [Wee 0; Wei 0];   % Matrix of E-to-E synaptic strengths
    Wei_mat = [0 Wie; 0 Wii];   % Matrix of E-to-I synaptic strengths
    
    % Set the stimuli such that stim_amp is the mean and the difference
    % (stim_amp1 - stim_amp2) is obtained from stim_diff
    stim_diff = 2*stim_amp_vec(parameter);           % Amplitude difference between 2 stimuli
    % stim_amp1 is amplitude of stimulus to Unit 1
    stim_amp1 = stim_amp_vec(parameter) + 0.5*stim_diff;
    % stim_amp2 is amplitude of stimulus to Unit 2
    stim_amp2 = stim_amp_vec(parameter) - 0.5*stim_diff;
    
    stim = zeros(N_Units, Nt);        % Initialize vector for the stimulus
    % then add stimulus to appropriate time points.
    t_off = t_on + stim_length;         % Time for stimulus offset
    stim(1,round(t_on/dt):round(t_off/dt)) = stim_amp1;
    stim(2,round(t_on/dt):round(t_off/dt)) = stim_amp2;
    
    % Now add noise to the stimulus while it is present
    % Loop through all the separate noise bins
    for bin = 1:N_noise_bins;
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
            stim_noise_amp*randn(2,1)*ones(1,i_off+1-i_on)/sqrt(noise_time_scale);
    end
    
    % Below is the formula for the firing rate response function that will
    % be used in both the steady state curves and the simulation.
    % This is a function of "x" where "x" is the value of the
    % variable (or vector of values) sent to the function.
    r_of_Vss = @(x) (x - Vth)./(tau*delta_V*(1-exp(-a*(x-Vth)/sigma_V)));
    
    % To prevent division by zero when x = Vth, use the limit
    % value:
    r_of_Vss_limit = 1./(tau*delta_V*a/sigma_V);
    
    r = zeros(N_Units,Nt);        % Initialize vector for firing rates
    sE = zeros(N_Units,Nt);       % Initialize vector for synaptic inputs
    sI = zeros(N_Units,Nt);       % Initialize vector for synaptic inputs
    D = ones(N_Units,Nt);         % Initialize vector for depression variable
    F = ones(N_Units,Nt);         % Initialize vector for facilitation variable
    
    % Now produce a vector of stimulus-independent noise for each cell
    V_Noise_vec = V_Noise*randn(N_Units,Nt)*sqrt(tau/dt);
    
    %% Now carry simulate the same model as a function of time
    % In models with multiple variables, each variable is updated using
    % the values of the set of variables on the previous time step
    
    % A while loop is used to exit the trial once threshold is
    % reached
    for i =2:Nt
        
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
        Vss_alt = Vss+V_Noise_vec(:,i);
        for cell = 1:N_Units
            if ( Vss_alt(cell) == Vth )
                r(cell,i) = r(cell,i-1)*(1-dt/tau) + ...
                    dt/tau*r_of_Vss_limit;
            else
                r(cell,i) = r(cell,i-1)*(1-dt/tau) + ...
                    dt/tau*r_of_Vss(Vss_alt(cell));
            end
        end
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
    
    %% Now plot the results
    if ( plot_all )
        figure(parameter)       % New figure for each stimulus
        clf
        % produce a solid gray rectangle across the region of the stimulus
        %                     fill([t_on t_on t_off t_off],[0 rmax rmax 0],[0.75, 0.75, 0.75], ...
        %                         'EdgeColor','none')
        hold on;        % do not overwrite the gray rectangle
        
        plot(tvec,r(1,:),'k');  % plot firing rate as a function of time
        plot(tvec,r(2,:),'k:');  % plot firing rate as a function of time
        xlabel('Time (sec)')        % Label x-axis on last trial
        ylabel('Rate (Hz)')
        axis([0.5 0.7 0 170])
    end
    %% Calculate oscillation frequency
    t_transient = 0.5;
    N_transient = round(t_transient/dt);
    r_calc1 = r(1,N_transient:end);
    r_calc2 = r(2,N_transient:end);
    N_calc = length(r_calc1);
    
    r_peak1 = max(r_calc1);
    r_min1 = min(r_calc1);
    r_peak2 = max(r_calc2);
    r_min2 = min(r_calc2);
    
    r_mid = (r_peak1+r_min1)/2;
    r_amp1(parameter) = r_peak1-r_min1;
    r_amp2(parameter) = r_peak2-r_min2;
    r_thresh_high = 0.75*r_peak1 + 0.25*r_min1;
    r_thresh_low = 0.25*r_peak1 + 0.75*r_min1;
    
    if ( r_calc1(1) > r_thresh_high )
        in_peak = 1;
    else
        in_peak = 0;
    end
    N_peaks = 0;
    t_cross = [];
    for i = 2:N_calc
        if ( ~in_peak && ( r_calc1(i) > r_thresh_high ) )
            in_peak = 1;
            N_peaks = N_peaks + 1;
            t_cross(N_peaks) = i;
        else
            if ( ( in_peak ) && ( r_calc1(i) < r_thresh_low ) )
                in_peak = 0;
            end
        end
    end
    Num_cycles(parameter) = N_peaks;
    if ( N_peaks > 1 )
        mean_f(parameter) = 1/(dt*(t_cross(N_peaks) - t_cross(1))/(N_peaks - 1));
    end
    mean_r1(parameter) = mean(r_calc1);
    mean_r2(parameter) = mean(r_calc2);
    
end;            % Go to next parameter set and produce a new trial

%% Now produce a power spectrum to find optimal frequency

f_vec = 0:0.5:100;          % Set of frequencies to test
Nf = length(f_vec);         % Number of distinct frequencies
for i = 1:Nf;               % Loop through set of frequencies
    f = f_vec(i);           % Value of frequency to test
    sinvec = sin(2*pi*f*tvec);      % sine function at that frequency
    cosvec = cos(2*pi*f*tvec);      % cosine function at that frequency
    
    % Now take the normalized dot produce of the sine function then the
    % cosine function with the firing rate vector. 
    sin_overlap = mean(sinvec.*r(1,:))/mean(r(1,:));
    cos_overlap = mean(cosvec.*r(1,:))/mean(r(1,:));
    
    % The power spectrum is given by the sum of the squares of the two
    % coefficients (sine coefficient and cosine coefficient). 
    % This essentially works because sin(x)*sin(x) + cos(x)*cos(x) = 1 for
    % all x.
    power(i) = (sin_overlap*sin_overlap + cos_overlap*cos_overlap);
end
figure(99)
plot(f_vec,power)

        

