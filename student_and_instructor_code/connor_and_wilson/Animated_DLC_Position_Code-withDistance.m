% Clear any existing variables and figures
clear variables;
close all;
camera_sFreq = 60;
smooth_win = 15;
% 4th and 10th xy in the csv

%Set WD
cd("Z:\NSB_2024\03_Mouse\Wilson\ConnorWiIson\Incoming\M521\Trial1\");

% Specify the filename of the CSV file
% filename = 'TRIMMED.M521_2024-07-04-R-Trial1_BaselineDLC_dlcrnetms5_SocialityJul6shuffle1_200000_el.csv';
filename = 'M521_2024-07-04-R-Trial1_FMRDLC_dlcrnetms5_SocialityJul6shuffle1_200000_el.csv';
TBL = readtable(filename,'NumHeaderLines',3);

% Extract data for mouse 1 and mouse 2
mouse1_x = TBL.x_3;
mouse1_y = TBL.y_3;
mouse2_x = TBL.x_9;
mouse2_y = TBL.y_9;

% Plot mouse 1 trajectory as a line in green
figure;
plot(mouse1_x, mouse1_y, 'g-', 'LineWidth', 2);
hold on;

% Plot mouse 2 trajectory as a line in purple
plot(mouse2_x, mouse2_y, 'k-', 'LineWidth', 2);

% Customize plot appearance
title('Trajectories of Mouse 1 and Mouse 2');
xlabel('X Axis');
ylabel('Y Axis');
legend('Mouse 1', 'Mouse 2');  % Add legend entries for mouse 1 and mouse 2



t_sec = TBL.coords/camera_sFreq;
clrs = lines(4);



%
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
plot(t_sec,D, 'b-', 'LineWidth', 2)


%average distance

% Calculate average distance between the two mice
avg_distance = mean(D(~BADROW));

% Plot average distance as a horizontal line across the plot
avg_distance_line = ones(size(t_sec)) * avg_distance;
plot(t_sec, avg_distance_line, 'r--', 'LineWidth', 2);








