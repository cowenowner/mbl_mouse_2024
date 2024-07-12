%%
clearvars; % Clears out variables so that we start fresh.
% Cowen code for loading.
cowen_code = 'C:\Users\Administrator\Documents\GitHub\mbl_mouse_2024\CowenLib';
matt_code = 'C:\Users\Administrator\Documents\GitHub\mbl_mouse_2024\vandermeerlab\code-matlab\shared';
addpath(genpath(cowen_code))

cd('D:\Data\M521_2024_07_05_Odor_Wilson_Connor_Trial_01_g0_a\M521_2024_07_05_Odor_Wilson_Connor_Trial_01_g0_imec0\kilosort_rawbin');
%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%
load AllSpikes;
% Load NIDQ wheel data.
[NIDQ] = NPXL_Extract_NIDQ(fullfile('..','..','M521_2024_07_05_Odor_Wilson_Connor_Trial_01_g0_t0.nidq.bin' ),1);
[NIDQ_fiber] = NPXL_Extract_NIDQ(fullfile('..','..','M521_2024_07_05_Odor_Wilson_Connor_Trial_01_g0_t0.nidq.bin' ),3);

% compute speed;
wheel_pos = NIDQ.data_V;

tmp_wheel_speed = [0 diff(wheel_pos)];
BADIX = abs(tmp_wheel_speed) > max(abs(tmp_wheel_speed))/2;
BADIX = movmean(BADIX,3)>0;
wheel_pos2 = wheel_pos;
wheel_pos2(BADIX) = nan;
wheel_pos2 = movmean(wheel_pos2,2000,'omitmissing');
wheel_speed = [0 diff(wheel_pos2)];

wheel_speed(BADIX) = nan;
wheel_speed = movmedian(wheel_speed,1500,'omitmissing');
wheel_speed = movmean(wheel_speed,2000,'omitmissing');
wheel_vel = wheel_speed ;
wheel_speed = abs(wheel_vel);
wheel_speed(wheel_speed>.00007) = nan;
wheel_speed = movmean(wheel_speed,2500,'omitmissing');

figure;plot(NIDQ.t_sec,wheel_speed)
yyaxis right
plot(NIDQ.t_sec,wheel_vel)

% Check to be sure the wheel data loaded...
figure; plot(NIDQ.t_sec, NIDQ.data_V(1,:)); title('wheel')
cd('D:\Data\M521_2024_07_05_Odor_Wilson_Connor_Trial_01_g0_a')
button_press_sec = load('M521_2024_07_05_Odor_Wilson_Connor_Trial_01_g0_tcat.nidq.xd_6_2_0.txt');
% make spikes ts -- should convert into function
nCells = length(SP);
end_time_sec = NIDQ.t_sec(end);

% sort by depth along probe
depths = [SP(:).neuropixels_depth_uM];
[~, sort_idx] = sort(depths, 'descend');
SP = SP(sort_idx);
% Plot the spiking activity along with the wheel postion.
T_sec = {SP.t_uS};
for ii = 1:length(T_sec)
    T_sec{ii} = T_sec{ii}/1e6;
end

[Q,~,Q_sec] = histcounts_cowen(T_sec,'binsize',.2);
% PCA in case you want to...
[PC,SC] = pca(Q);
SC1 = SC(:,1);
SC1 = movmean(SC1,30);
SC1_detrend = detrend(SC1');
% figure;plot(SC(:,1),SC(:,2),'.')
MUA = sum(Q,2);
MUA_smooth = movmean(MUA,20);
figure
ax(1) = subplot(2,1,1)
imagesc(Q_sec,[],Q');
plot_markers(button_press_sec)
xlabel('sec')
ax(2) = subplot(2,1,2);
plot(NIDQ.t_sec, NIDQ.data_V(1,:));
xlabel('sec')
ylabel('wheel V')
axis tight
plot_markers(button_press_sec)
linkaxes(ax,'x')

% Find the speeds and positions that match the spike matrix.
[SPD] = interp1(NIDQ.t_sec, wheel_speed, Q_sec);

figure
scatter(SPD,Q(:,1))
lsline

figure
scatter(SPD,SC1_detrend)
lsline

figure
subplot(1,2,1)
scatter(SPD,sum(Q,2))
lsline
subplot(1,2,2)
scatter(SPD,log10(sum(Q,2)))
lsline


figure
plot(NIDQ.t_sec, NIDQ.data_V(1,:));
xlabel('sec')
ylabel('wheel V')
axis tight
yyaxis right
plot(Q_sec,SC1_detrend');
hold on
plot(Q_sec,detrend(SC1_detrend'),'b-');



plot_markers(button_press_sec)


figure
ax(1) = subplot(2,1,1);
plot_raster(T_sec)
hold on
plot_markers(button_press_sec)
ylabel('Neuron ID')

ax(2) = subplot(2,1,2);
plot(NIDQ.t_sec, NIDQ.data_V(1,:));
xlabel('sec')
ylabel('wheel V')
axis tight
plot_markers(button_press_sec)

linkaxes(ax,'x')



% Matt code for visualization MultiRaster. But sadly need to clear Cowen path to avoid
% overlap of functions
rmpath(genpath(cowen_code))
addpath(genpath(matt_code))

myS = ts;
for iC = length(SP):-1:1

    SP(iC).t = SP(iC).t_uS * 10^-6;
    S.t{iC} = SP(iC).t;
    S.label{iC} = SP(iC).cluster_id;

end

%%
cfg = [];
evt = ts;
evt.t{1} = button_press_sec;
evt.label{1} = 'button';
MultiRaster(cfg, S)