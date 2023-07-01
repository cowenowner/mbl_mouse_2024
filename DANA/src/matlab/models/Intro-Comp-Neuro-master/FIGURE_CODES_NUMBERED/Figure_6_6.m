% Figure_6_6.m
% Simple firing-rate model of bistability.
% Uses sigmoidal f-I curve so that firing rate saturates.
% Synaptic feedback includes either depression or facilitation or a
% combination of the two.
% Plots firing rate and feedback curves to show fixed points then simulates
% the model with inputs.
%
% Trial 1 has depression only (facilitation factor set to 0).
% Trial 2 has facilitation only (depression time constant set to 0).
% Trial 3 has both depression and facilitation.
% This code was used to produce Figure 6.6 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%% Parameters for sigmoidal firing-rate response curve
rmax = 100;                 % Maximum rate
x_sigma = 0.1;              % Inverse of steepness of response
x_mid = 0.5;                % Default input for 1/2-max response (can vary)
tau = 0.010;                % Time constant for firing-rate to change

%% Parameters that vary across conditions are given in vectors, with each 
%  element of the vector to be used on one corresponding trial.
%  First trial is depression only.
%  Second trial is facilitation only.
%  Third trial has both facilitation and depression.

% Threshold for 1/2-activation is reduced in trial 2, increasing the firing 
% rate of neurons in the absence of input. This is possible as strong 
% facilitation allows for a weak baseline feedback to be used that ensures 
% spontaneous activity remains stable even at relatively high rates of 10Hz 
% or more. 
x_mid_vec = [0.45 0.25 0.45];   


%% Synaptic response parameters differ across the three trials

alpha_vec = [0.5 0.5 0.5];      % Baseline postsynaptic amplitude.
Wee_vec = [85 7 90];            % Strength of synaptic feedback.
fac_vec = [0 0.1 0.2];          % Facilitation factor (0 means no facilitation)
p0_vec = [1 0.1 0.1];           % Baseline vesicle release probability.

% Time constant for synaptic depression, zero means no depression
taud_vec = [0.15 0 0.15];       

tau_s = 0.002;                  % Synaptic time constant for AMPA receptors
tau_f = 0.5;                    % Faciliation time constant

%% Now set up vectors to generate the firing rate response to a range of 
%  inputs and the dependence of synaptic feedback on a range of rates.
S_range = [0:0.001:1];          % Range of synaptic inputs to use    
r_range = [0:0.01:rmax];        % Range of presynaptic rates to use.


%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf
%% Set up the simulation parameters 
dt = 0.0005;                % Time step to use
tmax = 2;                   % Maximum time for simulation
tvec = dt:dt:tmax;          % Vector of time points
Nt = length(tvec);          % Number of time points
Ntrials = 3;                % Number of different simulation trials

r = zeros(Ntrials,Nt);      % Initialize vector for firing rates
s = zeros(Ntrials,Nt);      % Initialize vector for synaptic inputs
D = ones(Ntrials,Nt);      % Initialize vector for depression variable
F = ones(Ntrials,Nt);      % Initialize vector for facilitation variable

t_on = 1;                   % Time for stimulus onset
t_off = 1.1;                % Time for stimulus offset
stim_amp = 0.2;             % Amplitude of stimulus
stim = zeros(1, Nt);        % Initialize vector for the stimulus
% then add stimulus to appropriate time points.
stim(round(t_on/dt):round(t_off/dt)) = stim_amp;    

% Loop through trials with different parameters on each trial
for trial = 1:Ntrials;         
    alpha0 = alpha_vec(trial);  % Value of baseline postsynaptic amplitude
    Wee = Wee_vec(trial);       % Value of total synaptic strength
    f_f = fac_vec(trial);       % Value of facilitation factor
    p0 = p0_vec(trial);         % Value of baseline release probability
    tau_d = taud_vec(trial);    % Value of time constant for depression
    x_mid = x_mid_vec(trial);   % Value to shift threshold of f-I curve
    
    % Below is the formula for the firing rate response function that will 
    % be used in both the steady state curves and the simulation.
    % This is the sigmoid function of "x" where "x" is the value of the
    % variable (or vector of values) sent to the function.
    r_of_s = @(x) rmax./(1 + exp(-(x-x_mid)/x_sigma));

    % Steady state formula for facilitation variable, F, depends on r.
    % See Appendix A and Chapter 5.
    fac_of_r = 1+(1/p0-1)*f_f*tau_f*r_range./(1+f_f*tau_f*r_range);
    
    % Below is release probability defined as pr = p0.F.
    pr = p0*fac_of_r;
    
    % Steady state formula for depression variable, D, depends on r.
    % See Appendix A and Chapter 5.
    dep_of_r = 1./(1+tau_d*pr.*r_range);
    
    % Postynaptic amplitude is defined as alpha = alpha_0*D
    alpha = alpha0*dep_of_r;
    
    % Steady state formula for the synaptic response fraction depends on r.
    % See Appendix A and Chapter 5.
    s_of_r = alpha.*pr.*r_range*tau_s./(1+alpha.*pr.*r_range*tau_s);
    
    % Total synaptic input is the synaptic response fraction scaled by
    % total synaptic strength, Wee via S(r) = Wee.s(r).
    S_in = Wee*s_of_r;
    
    % Generate firing rate as a function of total synaptic input
    r_out = r_of_s(S_range);
    
    % Set up the figure to plot S(r) and r(S).
    figure(1)    
    subplot('Position',[0.12 1.04-0.31*trial 0.33 0.22]) % Position of panel
    plot(S_range,r_out,'k');        % Plot r(S) (S varied r from formula)
    hold on;                        % Do not overwrite
    plot(S_in,r_range,'k--');       % Plot S(r) (r varied S from formula)
    axis([0 1 0 rmax])               % Set the range of each axis

    if ( trial == Ntrials )          % If it is not the last trial
        xlabel('Synaptic Input')    % Otherwise label the x-axis
    end
    ylabel('Rate (Hz)')             % Label the y-axis
    if ( trial == 2 )
        legend('rate', 'feedback')
    end

    %% Now carry simulate the same model as a function of time
    % In models with multiple variables, each variable is updated using
    % the values of the set of variables on the previous time step
    for i = 2:Nt;                               % Loop through time points

        % First update the firing rate using the equation
        % dr/dt = -(r + r(S + stim) )/tau 
        % where S = Wee.s and the value of s on the prior time-step is used
        r(trial,i) = r(trial,i-1)*(1-dt/tau) + ...
            dt/tau*r_of_s(Wee*s(trial,i-1)+stim(i));
        
        % Next update the facilitation variable which depends on firing
        % rate alone, using the equation
        % dF/dt = (1-F)/tau_F + f_f*(1/p0 - F)*r
        % Note that if f_f is zero then F is fixed at 1.
        F(trial,i) = F(trial,i-1) +  dt/tau_f*(1-F(trial,i-1)) + ...
            f_f*(1/p0-F(trial,i-1))*dt*r(trial,i-1);

        % Next update the depression variable, which depends on both firing
        % rate and the facilitation variable, using the equation
        % dD/dt = (1-D)/tau_D - pr.D.r (where pr = p0.F)
        if ( tau_d > 0 )        % tau_d = 0 ensures D = 1 always
            D(trial,i) = D(trial,i-1) +(1-D(trial,i-1))*dt/tau_d - ...
                p0*F(trial,i-1)*D(trial,i-1)*r(trial,i-1)*dt;
        else
            D(trial,i) = 1;     % No impact of depression if D = 1.
        end
        
        % Finally update the synaptic response variable,s, which depends on
        % all of rate, facilitation and depression variables, using the
        % equation
        % ds/dt = -s/tau_s + alpha.pr.r.(1-s)
        % where pr = p0.F and alpha = alpha_0*D
        s(trial,i) = s(trial,i-1)*(1-dt/tau_s) + ...
            dt*F(trial,i-1)*p0*D(trial,i-1)*alpha0*r(trial,i-1) ... 
            *(1-s(trial,i-1));
    end;        % go to next time point
    
    %% Now plot the results
    
    % Set up the position of the panel for plotting simulation results
    subplot('Position',[0.62 1.04-0.31*trial 0.36 0.22])
    
    % produce a solid gray rectangle across the region of the stimulus
    fill([t_on t_on t_off t_off],[0 rmax rmax 0],[0.75, 0.75, 0.75], ...
        'EdgeColor','none')
    hold on;        % do not overwrite the gray rectangle
    plot(tvec,r(trial,:),'k');  % plot firing rate as a function of time
    if ( trial == Ntrials )          % if not the last trial
        xlabel('Time (sec)')        % Label x-axis on last trial
    end
    ylabel('Rate (Hz)')
    axis([0 tmax 0 rmax])            % Set the range of each axis
    
end;            % Go to next parameter set and produce a new trial


% Finally label the subfigures
annotation('textbox',[0.00 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A1')
annotation('textbox',[0.5 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A2')
annotation('textbox',[0.00 0.68 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B1')
annotation('textbox',[0.5 0.68 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B2')
annotation('textbox',[0.00 0.38 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C1')
annotation('textbox',[0.5 0.38 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C2')
