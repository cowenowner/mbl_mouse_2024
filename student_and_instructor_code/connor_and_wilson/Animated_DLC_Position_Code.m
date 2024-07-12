%% Clear any existing variables and figures
clear variables;
close all;
camera_sFreq = 60;
smooth_win = 26;


%Set WD
cd("Z:\NSB_2024\03_Mouse\Wilson\ConnorWiIson\Incoming\M521\Trial1");

% Specify the filename of the CSV file
filename = 'M521_2024-07-04-R-Trial1_FMRDLC_dlcrnetms5_SocialityJul6shuffle1_200000_el.csv';
TBL = readtable(filename,'NumHeaderLines',3);
figure;
subplot(1,2,1)
plot(TBL.x_3,TBL.y_3,'b')
% hold on
% figure;
subplot(1,2,2)
plot(TBL.x_9,TBL.y_9,'r')



X{1} = TBL.x_3;
Y{1} = TBL.y_3;
X{2} = TBL.x_9;
Y{2} = TBL.y_9;
for ii = 1:2
    X{ii} = movmedian(X{ii},smooth_win,'omitmissing');
    Y{ii} = movmedian(Y{ii},smooth_win,'omitmissing');
end
BADROW = isnan(X{1}) | isnan(X{2}) | isnan(Y{1}) | isnan(Y{2});
% we could delete these rows. TBD

dX = X{1} - X{2};
dY = Y{1} - Y{2};
D = sqrt(dX.^2 + dY.^2);  % euclidian distance
groups= {'WT' 'FMR1'};

% x = TBL.x_3;
% y = TBL.y_3;
t_sec = TBL.coords/camera_sFreq;
clrs = lines(4);

figure
for ii = 1:2
    subplot(3,2,1)
    plot(X{ii},Y{ii},'.')
    hold on
    subplot(3,2,3:4)
    plot(t_sec,X{ii},'.','Color',clrs(ii,:))
    hold on
    plot(t_sec,Y{ii},'.','Color',clrs(ii,:))
end
subplot(3,2,5:6)
plot(t_sec,D)

subplot(3,2,1)
legend(groups)
%% new code for event detection%%%%%%%%%%%%%%%%%%%



% Calculate time associated with each frame
num_frames = height(TBL); % Total number of frames
time = (0:num_frames-1) / camera_sFreq; % Time in seconds for each frame

% Add 'Time' column to TBL
TBL.Time = time';



% Set interaction threshold (adjust as needed)
interaction_threshold = 25; % Adjust this threshold as needed

% Find indices where interaction occurs
interaction_indices = find(D < interaction_threshold);

% Initialize arrays to store start and end times
start_times = [];
end_times = [];

% Calculate 10-second windows centered on interaction periods
for i = 1:length(interaction_indices)
    % Calculate middle time of interaction
    middle_time = time(interaction_indices(i));
    
    % Calculate window start and end times
    window_start = middle_time - 5; % 5 seconds before middle_time
    window_end = middle_time + 4;   % 4 seconds after middle_time
    
    % Adjust window times to fit within session duration
    if window_start < 0
        window_start = 0;
    end
    if window_end > time(end)
        window_end = time(end);
    end
    
    % Store start and end times
    start_times = [start_times; window_start];
    end_times = [end_times; window_end];
end

% Create a table from start_times and end_times
results_table = table(start_times, end_times);

% Display the new table
disp('Non-overlapping Interaction Windows Table:');
disp(results_table);
% Sort table by start_times
results_table = sortrows(results_table, 'start_times');

%% Stephen's dumb verion
INTERACT_IX = D < 25;
smooth_win = camera_sFreq*.5;
INTERACT_IX = movmean(INTERACT_IX,smooth_win)>0;
df = [diff(INTERACT_IX); 0];
df(1) = 0;

start_IX = df  == 1;
end_IX = df == -1;
start_times_sec = t_sec(start_IX);
end_times_sec = t_sec(end_IX);
dur_interact_sec = end_times_sec - start_times_sec;

figure;histogram(dur_interact_sec); xlabel('Interaction time (sec)')





figure
plot(t_sec,D)
hold on
plot(t_sec,INTERACT_IX*100)

plot(results_table.start_times, zeros(size(results_table.start_times)),'g>')
plot(results_table.end_times, ones(size(results_table.end_times)),'r<')



% Initialize variables for merged interaction windows
merged_start = results_table.start_times(1);
merged_end = results_table.end_times(1);

% Ensure the merged window is 10 seconds long
if merged_end - merged_start < 10
    merged_end = merged_start + 10; % Extend the window to 10 seconds if it's shorter
end

merged_table = [];

% Merge overlapping or adjacent windows
for i = 2:height(results_table)
    if results_table.start_times(i) <= merged_end
        % Overlapping or adjacent, extend the merged window
        merged_end = max(merged_end, results_table.end_times(i));
        % Ensure the merged window does not exceed 10 seconds
        if merged_end - merged_start > 10
            merged_start = merged_end - 10;
        end
    else
        % Store the current merged window
        merged_table = [merged_table; merged_start, merged_end];
        % Start a new merged window
        merged_start = results_table.start_times(i);
        merged_end = results_table.end_times(i);
        % Ensure the new window is 10 seconds long if it's shorter
        if merged_end - merged_start < 10
            merged_end = merged_start + 10;
        end
    end
end

% Store the last merged window
merged_table = [merged_table; merged_start, merged_end];


% Create a new table with non-overlapping interaction periods
merged_results_table = array2table(merged_table, 'VariableNames', {'StartTime', 'EndTime'});

% Display the new table
disp('Non-overlapping Interaction Windows Table:');
disp(merged_results_table);

% Example: Save results to a CSV file
writetable(merged_results_table, 'WT_FMR1_non_overlapping_results_table.csv');
disp('Results saved to "non_overlapping_results_table.csv"');