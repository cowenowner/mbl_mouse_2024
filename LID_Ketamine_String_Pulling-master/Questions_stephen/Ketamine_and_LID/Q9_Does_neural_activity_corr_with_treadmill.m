function OUT = Q9_Does_neural_activity_corr_with_treadmill()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% this was tested on rat 320 day 3. C:\Users\Stephen Cowen\Box\Cowen Laboratory\Data\LID_Ketamine_Single_Unit_R56\Rat320\03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_IT = true;
SES = LK_Session_Info();
OUT = [];
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();
TRD = LK_Treadmill_Speed(EVT, META);
POS = LK_Load_and_Clean_POS;
IMU = LK_Load_and_Process_IMU;
% Load some LFP data...
load('./Processed_Data/best_channels_beta.mat','OUT')
%% Load the best non-reref.
if exist(fullfile('LFP',OUT.best_non_reref),'file')
    % load if stored locally.
    LFP = LK_Load_and_Clean_LFP('./LFP',OUT.best_non_reref);
else
    global DIRS
    lfp_dir = fullfile(DIRS.LFP_Dir,SES.rat_str,SES.session_str,'LFP');
    LFP = LK_Load_and_Clean_LFP(lfp_dir, OUT.best_non_reref);
end
filts = SPEC_create_filters({'beta'}, LFP.sFreq);
LFPfilt = filtfilt(filts{1},LFP.LFP);
% Let's determine if the IMU data correlates in some way to the treadmill
% data...
imu_in_trd = interp1(IMU.t_uS,IMU.speed,TRD.t_uS);
figure
IX1 = TRD.tread_direction == 1;
IX2 = TRD.tread_direction == -1;
plot(TRD.speed(IX1),imu_in_trd(IX1),'.',TRD.speed(IX2),imu_in_trd(IX2),'.')
lsline
xlabel('tread speed');ylabel('IMU speed');
% Let's determine if the PC of the spike data has some relationship with
% treadmill speed....
TSr = Restrict(TS,TRD.tread_intervals_uS);
bin_ms = 100;
[Q,Qt_uS] = Bin_ts_array(TSr,bin_ms*1000);
bin_t_uS = mean(Qt_uS,2);
Qf = conv_filter(Q,hanning(5));
[~,sc,lat] = pca(Qf);

scintrd = interp1_matrix(bin_t_uS,sc(:,1:2),TRD.t_uS);

figure
subplot(3,2,1:2)
imagesc(bin_t_uS,[],Qf')
subplot(3,2,3:4)
plot(bin_t_uS,sc(:,1),bin_t_uS,sc(:,2))
axis tight
subplot(3,2,5)
plot(TRD.speed(IX1),scintrd(IX1,1),'.',TRD.speed(IX2),scintrd(IX2,2),'.')
lsline
xlabel('tread speed');ylabel('sc spikes');

