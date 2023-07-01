% PSTH_osc.m
% PSTH_osc.m generates spikes as a random Poisson process with probability
% that oscillates sinusoidally with time. The random spike generation
% process is repeated 10 times to represent 10 trials of spike trains
% produced by a neuron. The spike trains are analyzed by multiple methods
% to produce different versions of the peri-stimulus time histogram (PSTH).
%
% Method 1: accumulation of spikes in successive time windows of 50ms.
% Method 2: accumulation of spikes in successive time windows of 200ms.
% Method 3: accumulation of spikes in a sliding window of 200ms.
% Method 4: accumulation of spikes filtered by a Gaussian function
%
% This code is used to produce Figure 3.2 in the textbook:
% An Introductory Course in Computational Neuroscience 
% by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;                          % Clear all variables from memory
tmax = 3;                       % Duration of each trial in sec
Ntrials = 10;                   % Number of trials
dt = 0.001;                     % Time-step in sec
t = 0:dt:tmax;                  % Vector of time points
Nt = length(t);                 % Number of time points

%% Set up the underlying firing rate of the inhomogeneous Poisson process
f = 1;                          % Frequency of oscillation in Hz
omega = 2*pi*f;                 % Convert to angular frequency
rdev = 6;                       % Amplitude of oscillation of firing rate
rmean = 12;                     % Mean firing rate

rate = rmean + rdev*sin(omega*t);   % Firing rate as a function of time

%% Now generate spikes according to that firing rate
spikes = zeros(Ntrials,Nt);     % Each row contains spikes for one trial

% For each trial generate a row vector that contains the number of spikes
% in each time bin. The probability of N spikes in time-bin (i) is given by
% the Poisson parameter, lambda = rate(i)*dt, as lambda^N*exp(-lambda)/N!
% While the time-bin is small enough that most bins contain 0 or 1 spike,
% it is possible that a few bins contain 2 or more spikes.
for trial = 1:Ntrials
    spikes(trial,:) = alt_poissrnd(rate*dt);
end

%% Method 1) Fixed bin-width of 50ms, consecutive bins
window_width = 0.05;                    % Set the window-width
Nwidth = round(window_width/dt);        % Number of time points in the window
Nwindows1 = round(Nt/Nwidth);            % Number of windows in the trial
tvals1 = dt*[Nwidth/2:Nwidth:Nt];       % Time points at centers of windows
psth1 = zeros(1,Nwindows1);              % Initialize vector of PSTH values
% Now loop through the windows, summing the spikes across all rows (:) and
% between all time-points in the window, such that the number of
% time-points in the sum is Nwidth, starting from the first time-point in
% the first window.
for i = 1:Nwindows1;
    psth1(i) = sum(sum(spikes(:,(i-1)*Nwidth+1:i*Nwidth)));
end
psth1 = psth1/(Ntrials*window_width);   % Convert to firing rate (spikes per trial per sec)

%% Method 2) Fixed bin-width of 200ms, consecutive bins
window_width = 0.2;                     % Set the window-width
Nwidth = round(window_width/dt);        % Number of time points in the window
Nwindows1b = round(Nt/Nwidth);          % Number of windows in the trial
tvals1b = dt*[Nwidth/2:Nwidth:Nt];      % Time points at centers of windows
psth1b = zeros(1,Nwindows1b);           % Initialize vector of PSTH values
% Now loop through the windows, summing the spikes across all rows (:) and
% between all time-points in the window, such that the number of
% time-points in the sum is Nwidth, starting from the first time-point in
% the first window.
for i = 1:Nwindows1b
    psth1b(i) = sum(sum(spikes(:,(i-1)*Nwidth+1:i*Nwidth)));
end
psth1b = psth1b/(Ntrials*window_width);   % Convert to firing rate (spikes per trial per sec)


%% Method 3) Sliding window of 200ms, all time-points
window_width = 0.2;                     % Set the window-width
Nwidth = round(window_width/dt);        % Number of time points in the window
Nwindows2 = Nt-Nwidth;                  % Number of windows in the trial
% In the sliding window, each time window just shifts by one time-step
% although all time-pints over the range of the window width are
% accumulated.
tvals2 = dt*[Nwidth/2+1:Nt-Nwidth/2];   % Time points at centers of windows
psth2 = zeros(1,Nwindows2);             % Initialize vector of PSTH values
% Now loop through the windows, summing the spikes across all rows (:) and
% between all time-points in the window, such that the number of
% time-points in the sum is Nwidth, starting from the first time-point in
% the first window.
for i = 1:Nwindows2
    psth2(i) = sum(sum(spikes(:,i:i+Nwidth-1)));
end
psth2 = psth2/(Ntrials*window_width);   % Convert to firing rate (spikes per trial per sec)


%% Method 4) Gaussian filter of 100ms width
window_width = 0.1;                     % smooth each spike by this amount
Nwidth = round(window_width/dt);        % No. of time-points in smoothing standard deviation
two_bw_sq = 2*window_width*window_width; % To be used in the Gaussian for every spike
tbins4 = dt*[0:Nt];                     % Time points to evaluate
Nbins4 = length(tbins4);                % No. of time points
psth4 = zeros(1,Nbins4);                % Initialize PSTH to zero
for trial = 1:Ntrials;                  % Loop through trials
    for i = find(spikes(trial,:))       % Time-bin numbers of the spikes
        % The prefactor spikes(trial(i)) could be greater than one if more
        % than one spike occurred in the time-bin of the spike. If
        % multiple spikes are not possible the factor can be removed.
        prefactor = spikes(trial,i);
        for bin = 1:Nbins4;             % For every time-point to plot
            tdiff = dt*(i-bin+0.5);     % Difference between time-point and spike-time
            % Now calculate the Gaussian function given the time-difference
            % between the time-point and spike-time divided by the
            % window-width. Each time-point accumulates contributions from
            % nearby spikes.
            psth4(bin) = psth4(bin) + prefactor*exp(-tdiff*tdiff/two_bw_sq);
        end
    end
end

% The following line normalizes the counts in each time-bin to produce a
% firing rate. The first line of normalization is the standard for a
% Gaussian with an extra division by Ntrials to get mean number of spikes
% per second per trial.
% The second line of normalization with the error function accounts for end
% effects. It evaluates to a factor of 1 (so has no impact) for time-bins
% far from the beginning and end of the trial. However, for time-bins near
% one end of the trial it can evaluate to 0.5 to account for the fact that,
% for example, the first time bin only receives contributions from spikes
% on one side,
for i = 1:Nbins4
    psth4(i) = psth4(i)/(Ntrials*sqrt(pi*two_bw_sq) ...
        *(0.5*(erf(tbins4(i)/window_width)+erf((tmax-tbins4(i))/window_width))));
end

%% Now set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf
subplot('Position',[0.18 0.84 0.8 0.14])    % First subfigure for spike trains
axis([0 tmax 0.5 Ntrials+0.5])    % x-axis is time (sec); y-axis is trial number

hold on;                % We want to draw all the spikes without erasing prior ones
box on;                 % Ensure a box is around the subfigure
% The next section will loop through trials plotting a small vertical line
% on a row corresponding to the trial number at x-values corresponding to
% the time of each spike.
for trial = 1:Ntrials
    for i = find(spikes(trial,:))
        % Plot spike-time on the x axis between two values on the y-axis
        % that correspond to the trial number.
        plot( [dt*i dt*i], [trial-0.4 trial+0.4],'k')
    end
end
set(gca,'YTick',[])             % Do not give y-axis values
ylabel({'Trials'; ''})          % Label y-axis (two lines used to align below)

figure(1)
hold on
subplot('Position',[0.18 0.65 0.8 0.14])    % Second subfigure

stairs(tvals1, psth1,'k','LineWidth',3)     % Plot PSTH as stair-case plot
axis([0 tmax 0 (rmean+rdev)*1.5])
ylabel('Rate (Hz)')

subplot('Position',[0.18 0.46 0.8 0.14])    % Third subfigure
stairs(tvals1b,psth1b,'k','LineWidth',3)      % Plot PSTH as stair-case plot
axis([0 tmax 0 (rmean+rdev)*1.5])
ylabel('Rate (Hz)')

subplot('Position',[0.18 0.27 0.8 0.14])    % Fourth subfigure
stairs(tvals2,psth2,'k','LineWidth',3)      % Plot PSTH as stair-case plot
axis([0 tmax 0 (rmean+rdev)*1.5])
ylabel('Rate (Hz)')

subplot('Position',[0.18 0.08 0.8 0.14])    % Fifth subfigure
stairs(tbins4,psth4,'k','LineWidth',3)      % Plot PSTH as stair-case plot
axis([0 tmax 0 (rmean+rdev)*1.5])
xlabel('Time (sec)')
ylabel('Rate (Hz)')

% Finally label each subfigure A-E
annotation('textbox',[0 0.98 0.03 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0 0.79 0.03 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0 0.60 0.03 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0 0.41 0.03 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','D')
annotation('textbox',[0 0.22 0.03 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','E')
