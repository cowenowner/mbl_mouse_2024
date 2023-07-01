% Tutorial_8_1.m
%
% This code produces a Hopfield-like network on a visible grid of units.
%
% This code was used to produce Figure 8.3 in the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
% The code can be run with different plasticity rules that predominantly
% differ in the means by which they produce synaptic depression.
% rule can take the values 1 to 5
% If rule = 1:
%   Potentiation only when both presynaptic and postsynaptic neurons fire
%   at a rate above the threshold.
%   To prevent runaway excitation at the end of each trial subtractive
%   normalization is carried out to ensure the mean of the weight matrix is
%   unchanged.
% If rule = 2:
%   Potentiation when both presynaptic and postsynaptic neurons fire at a
%   rate above the threshold.
%   Depression if either presynatpic or postsynaptic neuron fire at a rate
%   above the threshold while the other neuron fires at a rate below the
%   threshold.
% If rule = 3:
%   Potentiation when both presynaptic and postsynaptic neurons fire at a
%   rate above the threshold.
%   Depression if the presynatpic neuron fires at a rate above the
%   threshold while the postsynaptic neuron fires at a rate below the
%   threshold.
% If rule = 4:
%   Potentiation when both presynaptic and postsynaptic neurons fire at a
%   rate above the threshold.
%   Depression if the postsynatpic neuron fires at a rate above the
%   threshold while the presynaptic neuron fires at a rate below the
%   threshold.
% If rule = 5:
%   A quadratic dependence on postsynaptic firing rate and linear
%   dependence on presynaptic firing rate. Depression occurs for
%   postsynaptic rates below the threshold and potentiation for
%   postsynaptic rates above the threshold as in Figure 6.1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rule = 4;   % see preamble
rng(7)      % start the random number generator with a given seed

N1 = 17;    % size of edge of square array
viewmatrix = zeros(N1,N1); % square array for viewing
Ncells = numel(viewmatrix); % no. of cells in the array
Ntrials = 404;  % No. of trials for learning
inputstrength = 50; % input strength during training

Wmean = -0.3/Ncells;            % This mean is maintained with rule 1.
% Weight matrix W has W(i,j) as the connection strength
% from Unit i to Unit j.
W = zeros(Ncells,Ncells)+Wmean; % set initial connections to < 0.

%% Generate the patterns that will be used as the basis for a set of
%  stimuli. These stimulus patterns should not have substantial overlap.

Npatterns = 4; % Number of stimulus patterns to define and learn

% Each pattern initialized to zero at all entries and of the size of the
% matrix to view the data.
pattern1 = zeros(size(viewmatrix));
pattern2 = zeros(size(viewmatrix));
pattern3 = zeros(size(viewmatrix));
pattern4 = zeros(size(viewmatrix));

% see patterns by imagesc(pattern1) etc;
% First pattern is letter "X"
for i = 1:N1
    pattern1(i,i) = 1;                  % Top left to bottom right
    pattern1(i,N1+1-i) = 1;             % Top right to bottom left
end

%Second pattern is letter "Y"
for i = 1:(N1)/2-1;
    pattern2(i,i+2) = 1;                % Top left to center
    pattern2(i,i+3) = 1;                % Double the thickness
    pattern2(i,N1-1-i) = 1;             % Top right to center
    pattern2(i,N1-2-i) = 1;             % Double the thickness
end
for i = ceil((N1+2)/2)-3:N1
    pattern2(i,floor((N1+2)/2)) = 1;    % vertical line in center
end

% 3rd pattern is letter "Z"
for j = 1:N1
    pattern3(1,j) = 1;                  % Top horizontal line
    pattern3(N1,j) = 1;                 % Bottom horizontal line
    if ( j > 1 )
        pattern3(N1+2-j,j) = 1;          % Diagonal of the "Z"
    end
end

% 4th pattern is letter "O"
for i = 2:N1-1
    pattern4(i,2) = 1;                  % Top horizontal line
    pattern4(i,N1-1) = 1;               % Bottom horizontal line
end
for j = 2:N1-1
    pattern4(2,j) = 1;                  % Left vertical line
    pattern4(N1-1,j) = 1;               % Right vertical line
end

% Here we define a 3D array, to store all 4 of the 2D patterns stacked on
% top of each other.
all_patterns = zeros(N1,N1,Npatterns);
all_patterns(:,:,1) = pattern1;
all_patterns(:,:,2) = pattern2;
all_patterns(:,:,3) = pattern3;
all_patterns(:,:,4) = pattern4;

%% Simulation setup
dt = 0.001;         % time step for simulation
tau = 0.010;        % time constant for cells
tmax = 1;           % maximum time to wait
t = 0:dt:tmax;      %time vector
Nt = length(t);     % number of time points

%% Parameters for the sigmoidal firing-rate model
rmax = 50;              % maximum firing rate
gain = 1;               % steepness of sigmoid f-I curve
threshold = 10;         % input to achieve half maximum rate

%% Set parameters for plasticity rules and plasticity rates
switch rule
    case 1
        rate_thresh = rmax/2;       % threshold for potentiation
        epsilon = 0.1/Ncells;       % rate of plasticity
    case 2
        rate_thresh = rmax/2;       % threshold for potentiation
        epsilon = 0.1/Ncells;       % rate of plasticity
        LTD_factor = 0.15;          % relative rate for depression
    case 3
        rate_thresh = rmax/2;       % threshold for potentiation
        epsilon = 0.1/Ncells;       % rate of plasticity
        LTD_factor = 0.4;           % relative rate for depression
    case 4
        rate_thresh = rmax/2;       % threshold for potentiation
        epsilon = 0.1/Ncells;       % rate of plasticity
        LTD_factor = 0.25;          % relative rate for depression
    case 5
        rate_thresh = 0.9*rmax;     % threshold for potentiation
        epsilon = 0.005/Ncells;     % rate of plasticity
end

Wmax = 8/Ncells;     % Maximum connection strength
Wmin = -8/Ncells;    % Minimum connection strength

% Training is with degraded stimuli based on a particular stimulus pattern
% but with a probability for each individual input to be incorrect. The
% probability is given by the variable "errorprob_train" and in this case
% denotes the probability that a unit that should be "on" is not given
% input and that a unit that should be "off" is given input.
errorprob_train = 0.1;

for trial = 1:Ntrials                       % iterate through trials
    rate = zeros(Nt,Ncells);                % initialize rates to zero
    selectpattern = randi(Npatterns);   % choose a pattern randomly

    % Now set the chosen pattern to be the current trial's input pattern
    input_pattern = all_patterns(:,:,selectpattern);
        
    % The next line uses a vector of random numbers between 0 and 1, with
    % each random number corresponding to a unit. If the random number is
    % less than the error probability then that unit's input is flipped
    % from off to on (0 to 1) or on to off (1 to 0).
    % flipinputs is the set of units to have their inputs flipped
    flipinputs = find(rand(Ncells,1) < errorprob_train );
    % The following line switches 0s to 1s and 1s to 0s for the set of
    % units that have just been designated as errors.
    input_pattern(flipinputs) = 1-input_pattern(flipinputs);        
    
    % Now simulate through time with the degraded input pattern
    
    for i = 2:Nt                    % simulate network response to input
        if ( i < Nt / 2 )           % Only input for first 1/2 of trial
            % Input current to a unit is given by the external input (first
            % term) and the network feedback (second term)
            current = input_pattern(:)'*inputstrength + rate(i-1,:)*W;
        else;                           % Second half of the trial           
            current= rate(i-1,:)*W;     % Network feedback
        end
        % Sigmoidal f-I curve is used in the next line
        rss = rmax./(1+exp(-gain*(current-threshold))); 
        % A robust method for approaching input-dependent rate with time
        % constant "tau", i.e. using:
        % d(rate)/dt = (rss - rate)/tau
%        rate(i,:) = rss + (rate(i-1,:)-rss)*exp(-dt/tau); 
        rate(i,:) = rate(i-1,:) + (rss-rate(i-1,:))*dt/tau;
    end;        % Go to next time point
    viewmatrix(:) = rate(end,:);   % viewmatrix is purely for viewing
  
    %% Now use the firing rates across the trial to update the connection 
    %  strengths. Note that this is batch updating instead of continuous
    %  updating by altering W on each time point.
    
    switch rule
        case 1
            % + for presynaptic rate and postsynaptic rate above threshold
            dW = (double(rate'>rate_thresh))*(double(rate>rate_thresh));
        case 2
            % + for presynaptic rate and postsynaptic rate above threshold
            dW = (double(rate'>rate_thresh))*(double(rate>rate_thresh));
            % - for presynaptic rate above and postsynaptic rate below
            % threshold
            dW = dW - LTD_factor* ...
                (double(rate'>rate_thresh))*(double(rate<rate_thresh));
            % - for postsynaptic rate above and presynaptic rate below
            % threshold
            dW = dW - LTD_factor* ...
                (double(rate'<rate_thresh))*(double(rate>rate_thresh));
        case 3
            % + for presynaptic rate and postsynaptic rate above threshold
            dW = (double(rate'>rate_thresh))*(double(rate>rate_thresh));
            % - for presynaptic rate above and postsynaptic rate below
            % threshold
            dW = dW - LTD_factor* ...
                (double(rate'>rate_thresh))*(double(rate<rate_thresh));
        case 4;
            % + for presynaptic rate and postsynaptic rate above threshold
            dW = (double(rate'>rate_thresh))*(double(rate>rate_thresh));
            % - for postsynaptic rate above and presynaptic rate below
            % threshold
            dW = dW - LTD_factor* ...
                (double(rate'<rate_thresh))*(double(rate>rate_thresh));
        case 5
            % Quadratic dependence on postsynaptic rate, linear on
            % presynaptic rate
            dW = rate'*(rate.*(rate-rate_thresh)); % update for coactive cells
    end
    W = W+dW*epsilon*dt;    % Update connection strengths
    W = min(W,Wmax);        % maximum connection strength
    W = max(W,Wmin);        % minimum connection strength is inhibitory
    
    switch rule
        case 1
            % normalize connections to maintain the overall mean
            W = W - mean(mean(W))+Wmean;
    end
    
    % Plot progress every 100 trials
    if ( mod(trial,100) < 5)
        figure(selectpattern)
        subplot(2,1,1)
        imagesc(input_pattern)          % Input to network
        subplot(2,1,2)
        imagesc(viewmatrix)     % Response at end of trial
        drawnow
        caxis([0 rmax])
    end
    
end;        % Go to next trial

% Now move to testing memory of inputs in same manner, but no update of W
errorprob = 0.2; % produce degradation if non-zero
% Now test a pattern
inputstrength = inputstrength*0.5; % use weaker than original input
figure(99)
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
clf

for trial = 1:Npatterns
    rate = zeros(Nt,Ncells);
    selectpattern = trial;

    % Now set the chosen pattern to be the current trial's input pattern
    input_pattern = all_patterns(:,:,selectpattern);

    % Now corrupt the pattern by switching a fraction of the inputs
    flipinputs = find(rand(Ncells,1) < errorprob );
    input_pattern(flipinputs) = 1-input_pattern(flipinputs);
    
    figure(99)

    subplot('Position',[0.03+0.25*(trial-1) 0.52 0.2 0.35])
    imagesc(input_pattern);    % View input patterns
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    title(strcat(['Final Input ', ' ', num2str(selectpattern)]))
    colormap(gray)
    for i = 2:Nt
        if ( i < Nt / 2 ) % Only input for first 1/2 of trial
            current = input_pattern(:)'*inputstrength + rate(i-1,:)*W;
        else
            current= rate(i-1,:)*W;     % Network evolves due to internal structure
        end
        % Sigmoidal f-I curve is used in the next line
        rss = rmax./(1+exp(-gain*(current-threshold)));
        % Update with exponential Euler method
        rate(i,:) = rss + (rate(i-1,:)-rss)*exp(-dt/tau);
    end
    
    % viewmatrix is already defined as 2D the next line fills the elements
    viewmatrix(:) = rate(end,:);
     
    % Finally plot all data on one figure 
    figure(99)
    subplot('Position',[0.03+0.25*(trial-1) 0.04 0.2 0.35])
    imagesc(viewmatrix);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    title(strcat(['Final Response ', ' ', num2str(selectpattern)]))
    colormap(gray)
end