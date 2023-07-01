%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
GP = LK_Globals;
PLOT_IT = true;
mfile = 'Q2_Does_paw_track_corr_with_pos_imu_Ana';
% ana_dir
% ses_to_ana = 'Q2_Does_paw_track_corr_with_pos_imu';
ana_dir = 'C:\Users\Stephen Cowen\Dropbox\Foldershare\Analysis_Results_Dropbox\Q2_Does_paw_track_corr_with_pos_imu';
d = dir(fullfile(ana_dir,'Dset*.mat'));
% The following makes Matlab not do funning things with underscores in
% text.

ALL = [];
CORS = [];
sescnt = 1;
%%%%%%%%%%%%%%%%%%%%

for iF = 1:length(d)
    Dset = load(fullfile(ana_dir,d(iF).name));
    ALL(iF).rat = Dset.SES.rat;
    ALL(iF).session = Dset.SES.session;
    ALL(iF).mn_IMU_speed = Dset.mn_IMU_speed;
    ALL(iF).mn_POS_speed = Dset.mn_POS_speed;
    ALL(iF).mn_Rot_speed = Dset.mn_Rot_speed;
    ALL(iF).mn_Left_Paw_speed = Dset.mn_Left_Paw_speed;
    ALL(iF).mn_Nose_speed = Dset.mn_Nose_speed;
    ALL(iF).mn_Right_Paw_speed = Dset.mn_Right_Paw_speed;
    CORS(:,:,iF) = Dset.CORS;
end
% This converts the structre above to a useful table.
TBL = struct2table(ALL);
%% Summarize the statistics.

figure
subplot(4,1,1)
bar(TBL.session,TBL.mn_IMU_speed)
xlabel('Session')
ylabel('IMU Speed')

subplot(4,1,2)
bar(TBL.session,TBL.mn_POS_speed)
xlabel('Session')
ylabel('POS Speed')

subplot(4,1,3)
bar(TBL.session,TBL.mn_Rot_speed)
xlabel('Session')
ylabel('Rot Speed')

subplot(4,1,4)
plot(TBL.session,TBL.mn_Left_Paw_speed,'o-')
hold on
plot(TBL.session,TBL.mn_Right_Paw_speed,'*-')
plot(TBL.session,TBL.mn_Nose_speed,'p-')
legend('Left','Right','Nose'); legend boxoff
xlabel('Session')
ylabel('Paw Speed')

%% % Look at the correlations between all vbls...
C = nanmean(CORS,3);
figure
imagesc(C)
colorbar
set(gca,'XTick',1:Cols(CORS))
set(gca,'XTickLabel',Dset.CORS_lbls)
xtickangle(90)
set(gca,'YTick',1:Cols(CORS))
set(gca,'YTickLabel',Dset.CORS_lbls)