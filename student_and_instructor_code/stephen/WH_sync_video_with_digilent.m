% Synchronize the video with the mice.
% Cowen 2024
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% define the files.
digilent_file = 'C:\MBL DATA\test5\acq2.csv';
video_file = 'C:\MBL DATA\test5\recording_test5.mp4';
dlc_file = 'C:\MBL DATA\test1\recording_testdownsampledDLC_resnet50_hand-testJun30shuffle1_4000.csv';
USE_DIO = true;

% load the digilent data and get the timestamps for each video frame. 
DATA = readtable(digilent_file);
figure
plot(DATA.Time_s_, DATA.Channel1_V_, DATA.Time_s_,DATA.Channel2_V_)


if USE_DIO
    plot(DATA.Time_s_,DATA.DIO0)
    d = diff(DATA.DIO0);
else
    d = diff(DATA.Channel2_V_ > .5);
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
figure(101)
clf
subplot(2,2,3:4)
GIX = DATA.Time_s_ >= frame_time_sec(1) & DATA.Time_s_ <= frame_time_sec(end);
plot(DATA.Time_s_(GIX),DATA.Channel1_V_(GIX),DATA.Time_s_(GIX),DATA.Channel2_V_(GIX))
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
    % currAxes.Visible = "off";

    title(sprintf('Frame %d',frame_cnt))

    pause(1/v.FrameRate)

    subplot(2,2,2)
    % find the record closest to this frame.
    rec_ix = find(DATA.Time_s_ >= frame_sec,1,'first');
    plot(frame_sec,DATA.Channel1_V_(rec_ix),'b.')
    hold on
    plot(frame_sec,DATA.Channel2_V_(rec_ix),'r.')

    subplot(2,2,3:4)
    plot(frame_sec,DATA.Channel1_V_(rec_ix),'bo')
    hold on
    plot(frame_sec,DATA.Channel2_V_(rec_ix),'ro')
    
    frame_cnt = frame_cnt + 1;
end

%%
% Plot the data along with the video, frame by frame.

% If there is a DLC file, then load this and make a more fancy plot of the
% data with the movement.

dlc = readtable(dlc_file,'NumHeaderLines',2);
