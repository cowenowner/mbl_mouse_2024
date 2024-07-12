% Aligning photometry data for DA to events detected from Deep Lab Cut
% Author Connor Neifert and Stephen Cowen
% Date: 07/10/2024

clear variables;
close all;
camera_sFreq = 60; % Hz
smooth_win = 26;
bckWin = 2;     % time window to lookback before an event (seconds)
fwdWin = 10;     % time window to lookforward after an event (seconds)
interactionType = 'Trial 1 Midbody_Midbody_Fake';       % will also be used for output file naming
Dthreshold = 50;

data_dir = 'C:\Users\Administrator\Documents\Smelly\M521\DLC CSVs'; % directory for DLC .csv to be assessed
filename = 'M521_2024-07-04-R-Trial1_FakeDLC_dlcrnetms5_SocialityJul6shuffle1_200000_el.csv';    % filename of DLC .csv
dlc_fname = fullfile(data_dir, filename);   % assemble dlc filename
digilent_preprocessed_dir = 'C:\Users\Administrator\Documents\Smelly\M521\DigilentPreprocessed\'; % processed digilent data dir
digilent_preprocessed_fname = fullfile(digilent_preprocessed_dir,"M521_R_2024_07_04_Trial1_Fake_processed.mat");  % Processed digilent .mat fname
digilent_fname = ['C:\Users\Administrator\Documents\Smelly\M521\DigilentRaw\' ...
    'M521_R_2024_07_04_Trial1_Fake.csv']; % digilent raw data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the digilent data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DG_preproc = load(digilent_preprocessed_fname);
DG_data = readtable(digilent_fname);

frame_ix = find(diff(DG_data.DIO0*-1) > 0);
frame_times_sec = DG_data.Time_s_(frame_ix);
n_frames_loaded = length(frame_ix);

figure
plot(DG_data.Time_s_, DG_data.Channel1_V_)
hold on
plot(DG_data.Time_s_, DG_data.DIO0)
plot(frame_times_sec, zeros(size(frame_times_sec)),'r*')

%Set WD
cd(data_dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load DLC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TBL = readtable(dlc_fname,'NumHeaderLines',3);
frame_ix_in_DG = frame_ix(1:size(TBL,1));
frame_times_sec_in_DG = frame_times_sec(1:size(TBL,1));
TBL.DA = DG_data.Channel1_V_(frame_ix_in_DG);
TBL.DA_detrend = DG_preproc.detrend_60s(frame_ix_in_DG);
TBL.DA_z = DG_preproc.z_60s(frame_ix_in_DG);


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

% x = TBL.x;
% y = TBL.y;
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
yyaxis right
plot(t_sec,TBL.DA)

subplot(3,2,1)
legend(groups)
%
figure
plot(t_sec,D)
yyaxis right
plot(t_sec,TBL.DA)
% [pts, y] = ginput(6)

[r,p] = corr(D, TBL.DA, 'Rows', 'complete');
figure
scatter(D, TBL.DA)
lsline
xlabel('Distance to mouse')
ylabel('DA levels')

%% new code for event detection%%%%%%%%%%%%%%%%%%%
% Calculate time associated with each frame
num_frames = height(TBL); % Total number of frames
time = (0:num_frames-1) / camera_sFreq; % Time in seconds for each frame

% Add 'Time' column to TBL
TBL.Time = time';

%% Detect start and end times of interactions
INTERACT_IX = D < Dthreshold;
smooth_win = camera_sFreq*.5;
INTERACT_IX = movmean(INTERACT_IX,smooth_win)>0;
df = [diff(INTERACT_IX); 0];
df(1) = 0;

start_IX = df  == 1;
end_IX = df == -1;
start_times_sec = t_sec(start_IX);
end_times_sec = t_sec(end_IX);
dur_interact_sec = end_times_sec - start_times_sec;

while any(start_times_sec>=(t_sec(end)-fwdWin)) == 1;       % Remove any events with windows that will extend beyond the recording time
    start_times_sec(end) = [];
end
while any(start_times_sec<=(t_sec(1)-bckWin)) == 1;       % Remove any events with windows that will extend beyond the recording time
    start_times_sec = start_times_sec(2:end); 
end


figure;histogram(dur_interact_sec); xlabel('Interaction time (sec)')


%% visualize photometry data around each event initiation (start time)


colorShade = zeros(length(start_times_sec), 3);
colorShade(:,2) = linspace(1,0.5,length(start_times_sec));

for i = 1:length(start_times_sec)
    eventWin(i, :) = [start_times_sec(i)-bckWin*camera_sFreq, start_times_sec(i)+fwdWin*camera_sFreq];
    eventWin_IX(i,:) = [(find(t_sec == start_times_sec(i)))-bckWin*camera_sFreq, (find(t_sec == start_times_sec(i)))+fwdWin*camera_sFreq];
    
    eventWinDA(i,:) = TBL.DA_z(find(t_sec == start_times_sec(i)));

    tspan(i,:) = [t_sec(eventWin_IX(i,1):(eventWin_IX(i,2)))];      
    tspan_adjust = tspan(i,:)-tspan(i,1)-bckWin;                      % line up all traces so "event" begins at 0s
    DAspan(i,:) = [TBL.DA_detrend(eventWin_IX(i,1):(eventWin_IX(i,2)))];
    DAspan_z(i,:) = [TBL.DA_z(eventWin_IX(i,1):(eventWin_IX(i,2)))];

    DAspan_adjust(i,:) = DAspan(i,:)-DAspan(i,(bckWin*camera_sFreq+1));    % slide each trace so event starts at 0 fluorsence
    DAspan_z_adjust(i,:) = DAspan_z(i,:);% -DAspan_z(i,(bckWin*camera_sFreq+1));    
    DAdF_F_adjust(i,:) = (DAspan(i,:)-mean(DAspan(i:bckWin*camera_sFreq)))./mean(DAspan(i:bckWin*camera_sFreq));   %F/dF calculation
    
  
    figure(10)
    hold on
    plot(tspan_adjust, DAspan_adjust(i,:), 'Color',[colorShade(i,:)])        % plot the DA detrended trace for this iteration's interaction

    figure(11)
    hold on
    plot(tspan_adjust, DAspan_z_adjust(i,:), 'Color',[colorShade(i,:)])      % plot the DA z-score trace for this iteration's interaction


    figure(12)
    hold on
    plot(tspan_adjust, DAdF_F_adjust(i,:), 'Color',[colorShade(i,:)])        % plot the DA dF/F trace for this iteration's interaction

end

avgDAspan_adjust = mean(DAspan_adjust ,1);
avgDAspan_z_adjust = mean(DAspan_z_adjust ,1);
avgDAdF_F_adjust = mean(DAdF_F_adjust ,1);

% plotting extra deets for the graphs
figure(7)
hold on
plot(tspan_adjust,avgDAspan_adjust, 'k','LineWidth',2)
legend
xline(0)  % plot vertical line at x=0
c = "Dopamine_Detrended "+ interactionType + " Interactions";
title(c)
ylabel('Fluorsence (dV)'); xlabel('Time since event (s)')
hold off
% savefig(c)

figure(8)
hold on
plot(tspan_adjust,avgDAspan_z_adjust, 'k','LineWidth',2)
legend
xline(0)  % plot vertical line at x=0
c = "Dopamine z "+ interactionType + " Interactions";
title(c)
ylabel('Fluorsence (z-score)'); xlabel('Time since event (s)')
hold off
% savefig(c)


figure(9)
hold on
plot(tspan_adjust,avgDAdF_F_adjust, 'k','LineWidth',2)                     
legend
xline(0)   % plot vertical line at x=0
c = "Dopamine dF_F "+ interactionType + " Interactions";
title(c)
ylabel('dF/F (detrended)'); xlabel('Time since event (s)')
hold off
% savefig(c)