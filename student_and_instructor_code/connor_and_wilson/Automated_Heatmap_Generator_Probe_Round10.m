clear;
clc;

% Define folders and paths
signalFolderPath = 'C:\Users\hp\Desktop\PhD Stuff\MvdM Rotation\Data\SocialDA\Probe Test\Trimmed_Signals';
csvFolderPath = "C:\Users\hp\Desktop\PhD Stuff\MvdM Rotation\Data\SocialDA\Probe Test\Pos CSVs";
signal_heatMapFolderPath = "C:\Users\hp\Desktop\PhD Stuff\MvdM Rotation\Data\SocialDA\Probe Test\Probe Signal Heat Maps 10 Bin";
time_heatMapFolderPath = "C:\Users\hp\Desktop\PhD Stuff\MvdM Rotation\Data\SocialDA\Probe Test\Probe Time Heat Maps 10 Bin";

% Experimental Treatment schedule

% Get a list of .mat files in the folder
matFiles = dir(fullfile(signalFolderPath, '*.mat'));

% Loop through each .mat file
for fileIdx = 1:length(matFiles)
    % Load neural data struct
    neuralData = load(fullfile(signalFolderPath, matFiles(fileIdx).name));
    neuralData = neuralData.data;
    
    % Save file name as neuralfilename
    [~, fileName, ~] = fileparts(matFiles(fileIdx).name);
    neuralfilename = fileName;
    
    % Find matching CSV file
    csvFile = dir(fullfile(csvFolderPath, [neuralfilename '*DLC_resnet50_SocialDAMay7shuffle1_100000.csv']));
    
    % Load and process position data
    pos_data = csvread(fullfile(csvFolderPath, csvFile.name), 3, 0);

    % Remove rows with values < 0.85 in the 4th column
    pos_data(pos_data(:, 4) < 0.85, :) = [];

    pos_data(:, 5) = pos_data(:, 1) / 30; % Divide values in the 1st column by 30
    pos_data_times = pos_data(:, 5); % Save the 5th column as pos_data_times
    %pos_data(:, 5) = []; % Remove the 5th column temporarily
    
    % Initialize indices array
    indices = zeros(size(pos_data_times));
    
    % Find closest matches in neuralData.t to pos_data_times
    for i = 1:length(pos_data_times)
        [~, idx] = min(abs(neuralData.t - pos_data_times(i)));
        indices(i) = idx;
    end
    
    % Trim neuralData.t and neuralData.zdF
    neuralData.t_trimmed = neuralData.t(indices);
    neuralData.zdF_trimmed = neuralData.zdF(indices);
    
    % Determine the floor and ceiling values
    floorVal = neuralData.t_trimmed(1);
    ceilingVal = neuralData.t_trimmed(end);
    
    %Convert struct to matrix
    signal_data = neuralData.zdF_trimmed;
    signal_time = neuralData.t_trimmed;
    
    % Remove rows in pos_data based on floor and ceiling values
    %pos_data(pos_data_times < floorVal | pos_data_times > ceilingVal, :) = [];
    pos_data(pos_data_times < floorVal, :) = [];
    signal_data(pos_data_times < floorVal, :) = [];
    signal_time(pos_data_times < floorVal, :) = [];
    pos_data_times(pos_data_times < floorVal, :) = [];

    pos_data(pos_data_times > ceilingVal, :) = [];
    signal_data(pos_data_times > ceilingVal, :) = [];
    signal_time(pos_data_times > ceilingVal, :) = [];
    pos_data_times(pos_data_times > ceilingVal, :) = [];

    %Bin Position values to the nearest 10 number
    pos_data(:, 2:3) = round(pos_data(:, 2:3),-1);

    % Get the number of rows in the pos_data matrix
    numRows = size(pos_data, 1);
    
    % Create a column vector of ones with the same length as X
    onesColumn = ones(numRows, 1);
    
    % Add the ones column to the pos_data matrix
    pos_data = [pos_data, onesColumn];

    % Extract pos_data ones
    seconds = pos_data(:,6);
    
    % Create signal heat map
    heatMap = accumarray(pos_data(:, 2:3)/10, signal_data, [], @mean);

    % Create time heat map 
    heatMap_time = accumarray(pos_data(:, 2:3)/10, seconds, []);
    
    % Remove 0 values from the heat map

    heatMap(~heatMap) = NaN;
    heatMap_time(~heatMap_time) = NaN;

    % Replace NaN values with white color
    nanColor = 99999999; % Choose a unique placeholder value
    heatMap(isnan(heatMap)) = nanColor;
    heatMap_time(isnan(heatMap_time)) = nanColor;

    %Calculate limits
    parenthesesValue = str2double(extractBetween(fileName, '(', ')'));

    % Saving raw heatmaps
    outputFolder = 'C:\Users\hp\Desktop\PhD Stuff\MvdM Rotation\Data\SocialDA\Probe Test\Probe Signal Heat Maps 10 Bin';
    outputFolder_time = "C:\Users\hp\Desktop\PhD Stuff\MvdM Rotation\Data\SocialDA\Probe Test\Probe Time Heat Maps 10 Bin";

    % Create the full file paths
    heatMapFilePath = fullfile(outputFolder, [fileName '_heatmap.mat']);
    heatMapTimeFilePath = fullfile(outputFolder_time, [fileName '_heatmap_time.mat']);

    % Save the heatMap and heatMap_time as .mat files
    save(heatMapFilePath, 'heatMap');
    save(heatMapTimeFilePath, 'heatMap_time');

     %% Initialize the Mat_tx matrix
%     Mat_tx = [];
%     
%     % Check the value between parentheses and assign values to Mat_tx
%     if parenthesesValue == 383
%         Mat_tx = repmat(['GA', 'PR', 'WN'], 1, 4);
%     else
%         Mat_tx = repmat(['WN', 'GA', 'PR'], 1, 3);
%     end
% 
%     % Extract values 6-10 from fileIdx
%     dateSubstring = fileName(6:10);
%     
%     % Initialize the day_tx variable
%     day_tx = '';
%     
%     % Check values 6-10 of fileIdx and assign values to day_tx
%     if ismember(dateSubstring, {'05-02', '05-05', '05-08', '05-16'})
%         day_tx = Mat_tx(1:2);
%     elseif ismember(dateSubstring, {'05-03', '05-06', '05-09', '05-17'})
%         day_tx = Mat_tx(3:4);
%     else
%         day_tx = Mat_tx(5:6);
%     end
% 
%     % Initialize variables
%     xlow = 0;
%     xhigh = 0;
%     ylow = 0;
%     yhigh = 0;
%     
%     % Check the first value of day_tx and assign values accordingly
%     if strcmp(day_tx(1), 'G')
%         xlow = 35;
%         xhigh = 210;
%         ylow = 350;
%         yhigh = 560;
%     elseif strcmp(day_tx(1), 'P')
%         xlow = 150;
%         xhigh = 300;
%         ylow = 190;
%         yhigh = 350;
%     elseif strcmp(day_tx(1), 'W')
%         xlow = 220;
%         xhigh = 420;
%         ylow = 340;
%         yhigh = 560;
%     end

    %% Create signal figure
    signal_figure = figure;
    imagesc(heatMap);
    set(gca,'color',[0 0 0]);
    colormap([jet; 0 0 0]); % Add white color at the end of the colormap
    colorbar;
    caxis([-3 7]);
    xlim([4 42]);
    ylim([17 58]);
    xlabel('X Pixels');
    ylabel('Y Pixels');
    title('Heat Map');

    % Save the signal heat map as a JPEG file
    [~, heatMapFileName, ~] = fileparts(matFiles(fileIdx).name);
    signal_heatMapFilePath = fullfile(signal_heatMapFolderPath, [heatMapFileName '.jpg']);
    exportgraphics(gca,signal_heatMapFilePath,'Resolution',300);
    %saveas(signal_figure, signal_heatMapFilePath, 'jpg');
    %imwrite(heatMap, heatMapFilePath, 'JPEG');


    % Create time figure
    time_figure = figure;
    imagesc(heatMap_time);
    set(gca,'color',[0 0 0]);
    colormap([jet; 0 0 0]); % Add black color at the end of the colormap
    colorbar;
    caxis([0 600]);
    xlim([4 42]);
    ylim([17 58]);
    xlabel('X Pixels');
    ylabel('Y Pixels');
    title('Heat Map');
    
    % Save the time heat map as a JPEG file
    [~, heatMapFileName, ~] = fileparts(matFiles(fileIdx).name);
    time_heatMapFilePath = fullfile(time_heatMapFolderPath, [heatMapFileName '.jpg']);
    exportgraphics(gca,time_heatMapFilePath,'Resolution',300);
    %saveas(time_figure, time_heatMapFilePath, 'jpg');
    %imwrite(heatMap, heatMapFilePath, 'JPEG');
    
    %Print progress on loop
    fprintf('Computation and saving completed for: %s\n', fileIdx);

end
