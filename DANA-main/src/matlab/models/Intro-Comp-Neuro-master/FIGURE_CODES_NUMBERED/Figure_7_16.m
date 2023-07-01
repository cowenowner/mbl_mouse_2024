% Figure_7_16.m
% A code that produces a set of "avalanches" representing, for example,
% sequences of neural activation, with the number of neurons active in
% consecutive time-bins being the set of data. In this code an avalanche
% begins with a single neuron active. In subsequent time-bins, each active
% neuron can produce a number of newly active neurons, with the new number
% taken probabilistically from the Poisson distribution with mean of lambda
% that is by default equal to 1, to obtain the best approximation of
% criticality in this system.
%
% Active neurons do not themselves remain active in the next time-bin.
% When a time-bin occurs with no active neurons the avalanche has ended.
%
% In this code, we test the resulting distribution of avalanche sizes and
% durations for power-laws, and show the extent of scale-free behavior in
% the system.
%
% This code is used to produce Figure 7.16 in the book:
% "An Introductory Course in Computational Neuroscience"
% by Paul Miller, Brandeis University 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tic;                                        % Begin timing the code
% A huge number of trials is used to obtain good statistics. This number
% should be reduced if time for running the code is limited.
Ntrials = 250000;
max_duration = 5000;                    % Maximum avalanche duration
maximum_size = 5000;                    % Maximum instantaneous size

% The data of the instantaneous size of each avalanche at each point in
% time is stored in the array "Nfiring"
Nfiring = zeros(Ntrials,max_duration);

lambda = 1;                             % If not "1" should be very close

trial = 1;                              % Stores the current trial number
while ( trial <=  Ntrials);             % Proceed for Ntrials in total
    
    if ( mod(trial,200) == 0 )          % Just to show progress every 200 trials
        trial
    end
    
    time = 1;                           % First time-bin
    Nfiring(trial, time) = 1;           % One neuron fires in first time-bin
    
    % Continue till the maximum duration or no neurons fire in a time-bin
    while ( (time <= max_duration) && (Nfiring(trial,time) > 0 ) )
        
        current_size = Nfiring(trial,time);     % Number of neurons firing
        
        temp_sum = 0;                           % To accumulate number in next time-bin
        for n = 1:current_size;                 % Loop through all active neurons
            % Add the number newly activated by each active neuron
            temp_sum = temp_sum + poissrnd(lambda);
        end
        
        time = time + 1;                        % Go to next time-bin
        
        if ( temp_sum > maximum_size )          % If system size is reached
            time = max_duration+1;              % Done to end this trial
        else
            Nfiring(trial,time) = temp_sum;     % Record number of active neurons
        end
    end
    duration(trial) = time-1;                   % Duration of this trial
    % Only record avalanches of duration greater than 2, or complete by the
    % time the maximum duration is reached
    if ( ( time > 2 ) && ( time <= max_duration) )
        trial = trial + 1;                      % Go to next trial
    else
        Nfiring(trial,:) = 0;                   % Erase and repeat trial
    end
end

[new_durs, new_index ] = sort(duration);        % Sort into ascending order

Nfiring = Nfiring(new_index,:);                 % Sorted by duration
sum_size = sum(Nfiring,2);                      % Sorted by duration

%% The next section is used to find the mean avalanche size and profile at each duration
numdurs = 1;                                    % Number of distinct durations
dur = new_durs(1);                              % First distinct duration
temp_sum = 0;                                   % To find the mean "sum_size"
num_of_dur(1) = 1;                              % Initialize the number of avalanches of this duration
durations(1) = dur;                             % Record value of the first duration
for i = 1:Ntrials;                              % Loop through all trials
    if (new_durs(i) > dur )                     % If trial has a new duration 
        numdurs = numdurs + 1;                  % One more distinct duration
        temp_sum = sum_size(i);                 % Initialize temporary sum of sizes of new duration
        num_of_dur(numdurs) = 1;                % Initialize number of avalanches of new duration
        dur = new_durs(i);                      % Value of new duration
        durations(numdurs) = dur;               % Record value of new duration
    else
        temp_sum = temp_sum + sum_size(i);      % Accumulate the sum of total sizes of this duration
        num_of_dur(numdurs) = num_of_dur(numdurs) + 1;      % Number of avalanches of this duration
    end
    mean_sumsize(numdurs) = temp_sum/num_of_dur(numdurs);    % Calculate mean of prior duration
    
end

% In the next section, the "undersampling" of high durations is taken into
% account, so that where there is one avalanche over a range of durations,
% the frequency distribution is scaled as one divided by the range. 
% The range is considered to be from the midpoint of the previous value to
% the midpoint of the next value.
for i = 2:length(num_of_dur)-1
    num_of_dur(i) = 2*num_of_dur(i)/(durations(i+1)-durations(i-1));
end

[sorted_sizes, new_index ] = sort(sum_size);        % Sort into ascending order

%% The next section is used to find the number of avalanches of each distinct total size
numsizes = 1;                               % Number of distinct total sizes
this_size = sorted_sizes(1);                % First distinct total size
num_of_sizes(1) = 1;                        % No. of avalanches of this size
sizes(1) = this_size;                       % Current avalanche size
for i = 1:Ntrials;                          % Loop through all trials
    if (sorted_sizes(i) > this_size )       % If trial has a new total size
        numsizes = numsizes + 1;            % One more distinct total size
        num_of_sizes(numsizes) = 1;         % Initialize number of new total size
        this_size = sorted_sizes(i);        % Value of new total size
        sizes(numsizes) = this_size;        % Record value of new total size
    else
        num_of_sizes(numsizes) = num_of_sizes(numsizes) + 1;    % Number of avalanches of current size
    end
    
end

% In the next section, the "undersampling" of high total sizes is taken into
% account, so that where there is one avalanche over a range of sizes,
% the frequency distribution is scaled as one divided by the range. 
% The range is considered to be from the midpoint of the previous value to
% the midpoint of the next value.
for i = 2:length(num_of_sizes)-1
    num_of_sizes(i) = 2*num_of_sizes(i)/(sizes(i+1)-sizes(i-1));
end

%% Now calculate the exponents and check the scaling relation

% Fits to the logarithms of the form:
% log(y) = fo.p1 .log(x) + fo.p2 
% yield the coefficient y = A x^(fo.p1)
% so fo.p1 is the power-law coefficient


fit_dur = fit(log(durations(3:end)'),log(num_of_dur(3:end)'), ...
    'POLY1','Weight',max(num_of_dur(3:end),1));

fit_size = fit(log(sizes(3:end)'),log(num_of_sizes(3:end)'), ...
    'POLY1','Weight',max(num_of_sizes(3:end),1));

fit_size_dur = fit(log(durations(3:end)'),log(mean_sumsize(3:end)'),'POLY1', ...
    'Weight',num_of_dur(3:end)'); % Mean size vs duration

% The next line is a fit of size versus duration treating all avalanches
% separately rather than first taking the mean of each range. This method
% weights small avalanches (of which there are more of a given size) more
% strongly.

tau = -fit_size.p1                                % Power-law exponent "tau"
alpha = -fit_dur.p1                              % Power-law exponent "alpha"    
one_over_sigma_nu_z = fit_size_dur.p1               % Power-law exponent "1/(sigma.nu.z)"
one_over_sigma_nu_z_theory = (alpha-1)/(tau-1)  % Expected value of "1/(sigma.nu.z)"
gamma = one_over_sigma_nu_z-1;              % Value of gamma for scaling curves

%% In the next section the shapes of avalanches over time will be averaged 
%  together when they correspond to a particular range of durations.
% In order to average together avalanches with varying durations within a
% particular range, the duration of the avalanches are split into a number
% of bins, "Nplotbins".
 
Nplotbins = 20;         % Split each avalanche into this number of time-bins  

% Only consider avalanches with at least half as many time-bins as the
% number to be plotted (those that are too small are of indeterminate
% shape)
new_indices = find(new_durs >= Nplotbins/2);        % Trials to use
Ntrials_big = length(new_indices);                  % No. of trials to use
bigdurs = new_durs(new_indices);                    % Durations of these trials
bigNfiring = Nfiring(new_indices,:);                % Instantaneous sizes
scaled_bigNfiring = bigNfiring./(bigdurs.^gamma)';  % Now scale those sizes

unscaled_sizes = zeros(Ntrials_big,Nplotbins);      % To accumulate original sizes
scaled_sizes = zeros(Ntrials_big,Nplotbins);        % To accumulate scaled sizes

for i = 1:Ntrials_big;                              % Loop through trials to use
    
    % If avalanche has fewer time bins than number to be plotted
    if ( new_durs(i) <= Nplotbins )
        for j = 1:Nplotbins;            % For each bin to plot
            % Find the closest bin to use
            unscaled_sizes(i,j) = bigNfiring(i,ceil(bigdurs(i)*j/Nplotbins));
            scaled_sizes(i,j) = scaled_bigNfiring(i,ceil(bigdurs(i)*j/Nplotbins));
        end
    else;   % If avalanche has more time bins than the number to be plotted
        for j = 1:Nplotbins;            % For each bin to plot
            % Average over the corresponding range of time-bins in the
            % avalanche, matching the fractional duration
            unscaled_sizes(i,j) = mean(bigNfiring(i,1+ceil(bigdurs(i)*(j-1)/Nplotbins):ceil(bigdurs(i)*j/Nplotbins)) );
            scaled_sizes(i,j) = mean(scaled_bigNfiring(i,1+ceil(bigdurs(i)*(j-1)/Nplotbins):ceil(bigdurs(i)*j/Nplotbins)) );
        end
    end
    
end

numsizebins = 5;                % Number of segments of the data for separate averaging of the size-vs-time profiles
starts = [];                    % Starting trial for each segment
stops = [];                     % Last trial for each segment
binned_dur = [];                % Mean duration for each segment

for i = 1:numsizebins;          % Loop through each segment
    starts(i) = 1+round((i-1)*Ntrials_big/numsizebins);     % First trial of segment
    stops(i) = round(i*Ntrials_big/numsizebins);            % Last trial of segment
    binned_dur(i) = mean(bigdurs(starts(i):stops(i)));      % Mean duration of segment
end

unscaled_mean_sizes = zeros(numsizebins,Nplotbins);         % Initialize mean size vs time profiles
scaled_mean_sizes = zeros(numsizebins,Nplotbins);
for i = 1:numsizebins;                                      % Loop through segments of data
    unscaled_mean_sizes(i,:) = mean(unscaled_sizes(starts(i):stops(i),:));  % Mean of unscaled values
    scaled_mean_sizes(i,:) = mean(scaled_sizes(starts(i):stops(i),:));      % Mean of scaled values
end

% Now fit the data from the different segments assuming a power-law
% relationship between size and duration. This time find the exponent by
% fitting the means of each profile rather than from the previous exponent.
curve_fit = fit(log(binned_dur'),log(mean(unscaled_mean_sizes,2)),'POLY1');
gamma_opt = curve_fit.p1;           % Optimal value of gamma from fitting

% Now scale each avalanche by its duration using the optimal (fitted) value
% of gamma
optscaled_sizes = unscaled_sizes./( (bigdurs'.^gamma_opt)*ones(1,Nplotbins) );
for i = 1:numsizebins;          % Loop through all segments of data
    % Find the mean of the "optimally scaled" data for each segment
    optscaled_mean_sizes(i,:) = mean(optscaled_sizes(starts(i):stops(i),:));
end

%% Now plot the figure with 6 panels
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(11)
clf
% Panel A, log number of avalanches versus log duration
subplot('Position',[0.08 0.62 0.22 0.33])
plot(log(durations),log(num_of_dur),'.','Color',[0.6 0.6 0.6])
hold on
plot(log(durations),fit_dur.p2+fit_dur.p1*log(durations),'k')
set(gca,'XTick',[log(10) log(100) log(1000) ]);
set(gca,'XTickLabel',{'10^1' '10^2' '10^3' })
set(gca,'YTick',[log(1) log(100) log(10000)]);
set(gca,'YTickLabel',{'10^0' '10^2' '10^4'})
xlabel('Duration')
ylabel('No. of each value')
axis([log(2) log(5000) log(0.01) log(100000)])

% Panel B, log number of avalanches versus log total size
subplot('Position',[0.41 0.62 0.22 0.33])
plot(log(sizes),log(num_of_sizes),'.','Color',[0.6 0.6 0.6])
hold on
plot(log(sizes),fit_size.p2+fit_size.p1*log(sizes),'k')
set(gca,'XTick',[log(10) log(1000) log(100000) ]);
set(gca,'XTickLabel',{'10^1' '10^3' '10^5' })
set(gca,'YTick',[log(1) log(100)  log(10000)]);
set(gca,'YTickLabel',{'10^0' '10^2' '10^4' '10^5'})
axis([log(2) log(10000000) log(0.00001) log(100000)])
xlabel('Total Size')
ylabel('No. of each value')

% Panel C, log mean total size of avalanche versus duration
subplot('Position',[0.74 0.62 0.22 0.33])
plot(log(durations),log(mean_sumsize),'.','Color',[0.6 0.6 0.6])
hold on
plot(log(durations),fit_size_dur.p2+fit_size_dur.p1*log(durations),'k')
set(gca,'XTick',[log(10) log(100) log(1000) log(10000)]);
set(gca,'XTickLabel',{'10^1' '10^2' '10^3' '10^4'})
set(gca,'YTick',[log(10) log(1000) log(100000)]);
set(gca,'YTickLabel',{'10^1' '10^3' '10^5'})
xlabel('Duration')
ylabel('Mean Total Size')
axis([0 log(1e4) 0 log(1e7)])

% Create a set of colors so each segment of ordered data is a different
% shade of gray, the largest avalanches in the lightest gray.
% Note that a row [0 0 0] is black and [1 1 1] is white.
colorset = [0:numsizebins-1]'*ones(1,3)/numsizebins
set(groot,'defaultAxesColorOrder',colorset)

% Panel D. Plot the size versus time of the different segments of data
subplot('Position',[0.08 0.105 0.22 0.33])
plot([1:Nplotbins]'*binned_dur/Nplotbins, unscaled_mean_sizes')
xlabel('Time (bins)')
ylabel('Size')
axis([0 250 0 67])

% Panel E. Plot the scaled size versus normalized time for the different 
% segments of data, using the scaling obtained from Panel C
subplot('Position',[0.41 0.105 0.22 0.33])
plot([1:Nplotbins]'/Nplotbins,scaled_mean_sizes')
set(gca,'XTick',[0 1]);
set(gca,'YTick',[]);
xlabel('Time / Duration')
ylabel('Scaled size')

% Panel F. Plot the scaled size versus normalized time for the different 
% segments of data, using the scaling obtained from Panel C
subplot('Position',[0.74 0.105 0.22 0.33])
plot([1:Nplotbins]'/Nplotbins,opt2scaled_mean_sizes')
set(gca,'XTick',[0 1]);
xlabel('Time / Duration')
ylabel('Scaled size')
set(gca,'YTick',[]);

axis([0 1 0 0.5])

textinsert = strcat('\alpha = ',num2str(alpha, '%4.2f'));
annotation('textbox',[0.16 0.88 0.15 0.05],'String',textinsert,'LineStyle','none','FontSize',16,'FontWeight','bold')

textinsert = strcat('\tau = ',num2str(tau, '%4.2f'));
annotation('textbox',[0.5 0.88 0.15 0.05],'String',textinsert,'LineStyle','none','FontSize',16,'FontWeight','bold')

textinsert = strcat('1/\sigma\nuz = ',num2str(fit_size_dur.p1, '%4.2f'));
annotation('textbox',[0.74 0.88 0.25 0.05],'String',textinsert,'LineStyle','none','FontSize',16,'FontWeight','bold')

textinsert = strcat('\gamma = 1/\sigma\nuz - 1 = ',num2str(gamma, '%4.2f'));
annotation('textbox',[0.41 0.12 0.25 0.05],'String',textinsert,'LineStyle','none','FontSize',16,'FontWeight','bold')

textinsert = strcat('\gamma^{opt} = ',num2str(gamma_opt, '%4.2f'));
annotation('textbox',[0.79 0.135 0.25 0.05],'String',textinsert,'LineStyle','none','FontSize',16,'FontWeight','bold')

annotation('textbox',[0.00 0.96 0.05 0.05],'String','A','LineStyle','none','FontSize',18,'FontWeight','bold')
annotation('textbox',[0.33 0.96 0.05 0.05],'String','B','LineStyle','none','FontSize',18,'FontWeight','bold')
annotation('textbox',[0.66 0.96 0.05 0.05],'String','C','LineStyle','none','FontSize',18,'FontWeight','bold')
annotation('textbox',[0.00 0.43 0.05 0.05],'String','D','LineStyle','none','FontSize',18,'FontWeight','bold')
annotation('textbox',[0.33 0.43 0.05 0.05],'String','E','LineStyle','none','FontSize',18,'FontWeight','bold')
annotation('textbox',[0.66 0.43 0.05 0.05],'String','F','LineStyle','none','FontSize',18,'FontWeight','bold')

toc
save('avalanches_lambda_1.mat','-v7.3')


