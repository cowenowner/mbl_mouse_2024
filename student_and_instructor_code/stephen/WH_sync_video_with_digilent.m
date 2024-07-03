% Synchronize the video with the mice.
% Cowen 2024
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% define the files.
close all
digilent_file = 'D:\M520\2024-07-02\Trial 1\M520_L_2024_07_02_Trial1_Baseline.csv';
video_file = 'D:\M520\2024-07-02\Trial 1\M520_L_2024_07_02_Trial1_Baseline.mp4';
dlc_file = []; % The deep lab cut csv file if it exists
out_video_dir = "C:\Temp\";
vid_limits_x = [200 450];
vid_limits_y = [200 420];
% dlc_file = 'C:\MBL DATA\test1\recording_testdownsampledDLC_resnet50_hand-testJun30shuffle1_4000.csv';
USE_DIO = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,digilent_fname] = fileparts(digilent_file);
out_video_file = fullfile(out_video_dir,digilent_fname);
% load the digilent data and get the timestamps for each video frame. 
DATA = readtable(digilent_file);
DATA.smoothed = movmean(DATA.Channel1_V_,2000);

figure
plot(DATA.Time_s_, DATA.smoothed)


if USE_DIO
    plot(DATA.Time_s_,DATA.DIO0)
    d = diff(DATA.DIO0);
else
    d = diff(DATA.smoothed > .5);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% Validate and compute sample rates and num frames.
up_ix = find(d>0);
frame_time_sec = DATA.Time_s_(up_ix);
n_frames = length(frame_time_sec);
frame_interval_sec = median(diff(frame_time_sec));
video_sFreq = 1/frame_interval_sec;

fprintf('%d frames, %1.2f Hz\n',n_frames, video_sFreq)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% load the video file. Check to ensure the number of frames in the video
% data corresponds to the number of frames in the digilent file.
v = VideoReader(video_file);
vw = VideoWriter(out_video_file);
open(vw)

figure(101)
clf
subplot(2,2,3:4)
GIX = DATA.Time_s_ >= frame_time_sec(1) & DATA.Time_s_ <= frame_time_sec(end);
plot(DATA.Time_s_(GIX),DATA.smoothed(GIX))
axis tight
hold on

frame_cnt = 1;
while hasFrame(v)
    frame_sec = frame_time_sec(frame_cnt);

    subplot(2,2,1)
    % currAxes = axes;

    vidFrame = readFrame(v);
    % image(vidFrame,"Parent",currAxes)
    image(vidFrame)
    xlim(vid_limits_x)
    ylim(vid_limits_y)
    % currAxes.Visible = "off";

    title(sprintf('Frame %d Sec %1.4f',frame_cnt,frame_sec))

    % pause(1/v.FrameRate)

    subplot(2,2,2)
    % find the record closest to this frame.
    rec_ix = find(DATA.Time_s_ >= frame_sec,1,'first');
    plot(frame_sec,DATA.smoothed(rec_ix),'b.')
    hold on
    title(digilent_fname)

    subplot(2,2,3:4)
    plot(frame_sec,DATA.smoothed(rec_ix),'bo','MarkerSize',3)
    
    frame = getframe(gcf);
    writeVideo(vw,frame)

    frame_cnt = frame_cnt + 1;

    if mod(frame_cnt,200) == 0
        fprintf('%d ',frame_cnt)
    end

end
close(vw)
close(v)

%%
% Plot the data along with the video, frame by frame.

% If there is a DLC file, then load this and make a more fancy plot of the
% data with the movement.

dlc = readtable(dlc_file,'NumHeaderLines',2);
