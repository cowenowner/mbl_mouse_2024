function OUT = Q2_Does_paw_track_corr_with_pos_imu()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run this function in the data directory.
% Determine if paw tracking is correlating with IMU and POS
% This just analyzes behavior, not neural data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
global DIRS % Assumes this variable was set up already. Do this in LK_Analyze_sessions (or somewhere).
GP = LK_Globals;

OUT = [];
close all;
PLOT_IT = true;
SES = LK_Session_Info();
OUT.SES = SES;
Behavior_sFreq = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[MT, good_string_pull_intervals_uSec, PAW, ROT, POS, IMU] = LK_Combine_All_String_Pull_Motion_To_Table(Behavior_sFreq, true);

% Look at the correlation bewteen relevant movement parameters...
CORS_lbls = MT.Properties.VariableNames(2:end);
CORS = corr(Z_scores(MT{MT.Is_pulling,2:end}), 'type', 'Pearson','Rows','pairwise');

pull_times_sec = ROT.t_uSec(diff(ROT.Speed) > 30)'/1e6;
pull_times_sec = Restrict(pull_times_sec,good_string_pull_intervals_uSec/1e6);
neg_pull_times_sec = ROT.t_uSec(diff(ROT.Speed)< -30)/1e6;
neg_pull_times_sec = Restrict(neg_pull_times_sec,good_string_pull_intervals_uSec/1e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look for relationships between the behavioral data...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ROT.POS = interp1(POS.Time_uS,POS.COM_xy,ROT.t_uSec(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Send important things out for analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT.mn_Rot_speed = nanmean(MT.Rot_speed);
OUT.mn_IMU_speed = nanmean(MT.IMU_speed);
OUT.mn_POS_speed = nanmean(MT.POS_speed);
OUT.mn_Right_Paw_speed = nanmean(MT.Right_speed);
OUT.mn_Left_Paw_speed = nanmean(MT.Left_speed);
OUT.mn_Nose_speed = nanmean(MT.Nose_speed);
OUT.CORS = CORS;
OUT.CORS_lbls = CORS_lbls;

if PLOT_IT
    a = [min(ROT.t_uSec(ROT.GIX)/60e6) max(ROT.t_uSec(ROT.GIX)/60e6)];
    figure
    subplot(4,1,1)
    plot(ROT.t_uSec(ROT.GIX)/60e6,ROT.Speed(ROT.GIX));
    xlabel('min')
    hold on
    plot(pull_times_sec/60,ones(size(pull_times_sec))*nanmin(ROT.Speed(ROT.GIX)),'g>')
    % plot(neg_pull_times_sec/60,ones(size(neg_pull_times_sec))*nanmin(ROT.Speed(ROT.GIX)),'r<')
    
    title([SES.title_str ' Rotary Encoder Speed'])
    set(gca,'XLim',a)
    
    subplot(4,1,2)
    plot(POS.Time_uS(POS.GIX)/60e6,POS.Green_xy(POS.GIX,1),'g.')
    hold on
    plot(POS.Time_uS(POS.GIX)/60e6,POS.Green_xy(POS.GIX,2),'c.')
    title('XY Tracking Green')
    set(gca,'XLim',a)
    subplot(4,1,3)
    plot(IMU.t_uS(IMU.GIX)/60e6,IMU.data_V(IMU.GIX,1),'.')
    hold on
    plot(IMU.t_uS(IMU.GIX)/60e6,IMU.data_V(IMU.GIX,2),'.')
    title('IMU')
    set(gca,'XLim',a)
    if ~isempty(PAW)
        subplot(4,1,4)
        plot(PAW.Time_uSec(PAW.GIXr)/60e6,PAW.Right_x(PAW.GIXr),'.')
        hold on
        plot(PAW.Time_uSec(PAW.GIXr)/60e6,PAW.Right_y(PAW.GIXr),'.')
        plot(PAW.Time_uSec(PAW.GIXr)/60e6,PAW.Left_x(PAW.GIXr),'.')
        plot(PAW.Time_uSec(PAW.GIXr)/60e6,PAW.Left_y(PAW.GIXr),'.')
        ylabel('xy paw')
        title('Paw track')
        
        %         yyaxis right
        %         plot(PAW.Time_uSec(PAW.GIXr)/60e6,PAW.Right_x_to_nose(PAW.GIXr),'.')
        %         hold on
        %         plot(PAW.Time_uSec(PAW.GIXr)/60e6,PAW.Right_y_to_nose(PAW.GIXr),'.')
        %         ylabel('to NOSE')
        set(gca,'XLim',a)
        
        xlabel('min')
        
    end
    
    set(gcf,'Position',[ 474.6 33.8 560 836])
end

