function OUT = Q1_Do_neurons_correlate_with_string_pulling_ensemble()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run this function in the session directory.
% Determine if single-units oscillate more after ketamine injection.
% Compare autocorrs before and after injection.
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = [];
GP = LK_Globals;

Behavior_sFreq = 150;
neuron_quality_threshold = 2;
small_binsize_ms = 5;
params = 25:50:2000;
vbls = {'Right_speed' 'Left_speed' 'Nose_speed' 'Left_acc' 'Right_acc'  'Right_x_to_nose' 'Left_x_to_nose' 'Rot_speed' 'Rot_acc' 'IMU_speed' 'Right_y_d1' 'Left_y_d1'};

PLOT_IT = true;
SES = LK_Session_Info();
OUT.SES = SES;
OUT.aborted = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and create a magnificent table of motion data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MT, string_pull_intervals_uSec] = LK_Combine_All_String_Pull_Motion_To_Table(Behavior_sFreq, false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load spikes and restrict to the times of interest.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SP,TS_uS] = LK_Load_Spikes(neuron_quality_threshold,string_pull_intervals_uSec([1 end]));
OUT.n_neurons = length(SP);
if length(SP) < 4
    OUT.aborted = true;
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a reduced movement segment of just pulling times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MTp = MT(MT.Is_pulling,:);
% behavior_bin_size_ms = 1000/Behavior_sFreq;
cond = {'Rsq','Rsq_sh'};

for iCond = 1:2
    for iPar = 1:length(params)
        % smooth_win_ms = 200;
        smooth_win_ms = params(iPar);
        smooth_bins = round(smooth_win_ms/small_binsize_ms);
        switch cond{iCond}
            case 'Rsq'
                [NT,~,NT_t_uS] = Bin_and_smooth_ts_array(TS_uS, small_binsize_ms*1e3, smooth_bins);
            case 'Rsq_sh'
                [NT,~,NT_t_uS] = Bin_and_smooth_ts_array(TS_uS, small_binsize_ms*1e3, smooth_bins, true);
        end
        NT_MT = interp1(NT_t_uS,NT,MTp.t_uSec,'nearest');

        [~,NT_MT_pc] = pca(NT_MT);
        mxpc = min([Cols(NT_MT_pc) 5]);
        NT_MT_pc = NT_MT_pc(:,1:mxpc);
        for iY = 1:length(vbls)
            y = MTp.(vbls{iY});
            lm = fitlm(NT_MT_pc,y,'linear');% 'RobustOpts','on');
            OUT.(cond{iCond})(iY,iPar) = lm.Rsquared.Adjusted;
            OUT.([cond{iCond} '_p'])(iY,iPar) = lm.coefTest;
        end
    end
end
OUT.Rsq(OUT.Rsq_p > 0.05) = nan;
OUT.smooth_win_ms = params;
OUT.variables = vbls;
if PLOT_IT
    figure
    imagesc(OUT.smooth_win_ms,[],OUT.Rsq - OUT.Rsq_sh)
    set(gca,'YTick',1:length(OUT.variables))
    set(gca,'YTickLabel',OUT.variables)
    xlabel('msec')
    colorbar
end
