%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q13_Does_FR_correlate_with_movement_Ana_Abhi';

GP.Analysis_Dir = 'E:\Temp\TempAnaResults\Q13_Does_FR_correlate_with_movement_during_ketamine_Abhi_4.3.23';

PLOT_IT = false;
% adir = fullfile(GP.Analysis_Dir,ses_to_ana);
adir = fullfile(GP.Analysis_Dir);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end
cm = lines(5);

ses = {};

%%%%%%%%%%%%%%%%%%%%
NRN = [];
ALL_AC_base = [];
ALL_WV = [];
ALL_IMU = [];

nrn_cnt = 1;
IMU_cnt = 1;
for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));
    if categorical(Dset.SES.Drugs) == 'LDOPA&Ketamine'
        
        ALL_Rat(IMU_cnt).Rat = Dset.SES.Rat; 
        if Dset.SES.RatType == 'SHAM'
            ALL_Rat(IMU_cnt).Group = 'SHAM';
        elseif Dset.SES.RatType == '6ODHA_LID'
            ALL_Rat(IMU_cnt).Group = '6OHDA_LID';
        end
        
        ALL_IMU_bin(IMU_cnt,:) = Dset.IMU_All_bins_uS';        
        ALL_IMU(IMU_cnt,:) = Dset.IMU_Speed_All';
        
        IMU_cnt = IMU_cnt + 1;
    
    else
    end
end

        
%         ALL_IMU(IMU_cnt).Session = Dset.SES.Session;
%         ALL_IMU(IMU_cnt).Drugs = Dset.SES.Drugs;
%         ALL_IMU(IMU_cnt).IMU_All = Dset.IMU_Speed_All;
%         ALL_IMU(IMU_cnt).IMU_Base = Dset.IMU_Speed_base;
%         ALL_IMU(IMU_cnt).IMU_Post1 = Dset.IMU_Speed_post1;
%         ALL_IMU(IMU_cnt).IMU_Post2 = Dset.IMU_Speed_post2;
%         ALL_IMU(IMU_cnt).IMU_Post3 = Dset.IMU_Speed_post3;
%         ALL_IMU(IMU_cnt).IMU_All_bins_uS = Dset.IMU_All_bins_uS;
%         ALL_IMU(IMU_cnt).IMU_Base_bins_uS = Dset.IMU_base_bins_uS;
%         ALL_IMU(IMU_cnt).IMU_Post1_bins_uS = Dset.IMU_psot1_bins_uS;
%         ALL_IMU(IMU_cnt).IMU_Post2_bins_uS = Dset.IMU_post2_bins_uS;
%         ALL_IMU(IMU_cnt).IMU_Post3_bins_uS = Dset.IMU_post3_bins_uS;
%         ALL_IMU(IMU_cnt).binsize_ms = Dset.binsize_Q_S_ms;
        
%         for ii = 1:length(Dset.FRates_base)
%             
%             All_TS_al{nrn_cnt} = Dset.TSaligned{ii};
%             
%             Inj_uS = Dset.TM.EventStartEndUsec(1,2);
%             %             if categorical(Dset.SES.Drugs) == 'Saline&Ketamine'
%             %                 NRN(nrn_cnt).Inj_uS = Inj_uS;
%             %             else
%             Inj_uS = Inj_uS+(61*60e6);
%             NRN(nrn_cnt).Inj_uS = Inj_uS;
%             %             end
%             
%             NRN(nrn_cnt).NeuronID = Dset.SP(ii).fname;
%             NRN(nrn_cnt).Rat = Dset.SES.Rat;
%             NRN(nrn_cnt).Session = Dset.SES.Session;
%             if Dset.SES.RatType == 'SHAM'
%                 NRN(nrn_cnt).Group = 'SHAM';
%             elseif Dset.SES.RatType == '6ODHA_LID'
%                 NRN(nrn_cnt).Group = '6OHDA_LID';
%             end
%             NRN(nrn_cnt).Drugs = categorical(Dset.SES.Drugs);
%             if isempty(Dset.SP(ii).Depth_uM)
%                 dp = nan;
%             elseif Dset.SP(ii).Depth_uM ==0
%                 dp = nan;
%             else
%                 dp = Dset.SP(ii).Depth_uM;
%             end
%             NRN(nrn_cnt).Depth_uM = dp;
%             if dp < 2800
%                 NRN(nrn_cnt).BrainRegion = 'M1';
%             else
%                 NRN(nrn_cnt).BrainRegion = 'Striatum';
%                 NRN(nrn_cnt).M1Layer = 'Str';
%             end
%             
%             if dp < 500
%                 NRN(nrn_cnt).M1Layer = '123';
%             elseif dp >= 500 && dp < 900
%                 NRN(nrn_cnt).M1Layer = '4';
%             elseif dp >= 900 && dp < 2000
%                 NRN(nrn_cnt).M1Layer = '56';
%             elseif dp >= 2000 && dp < 2800
%                 NRN(nrn_cnt).M1Layer = 'M1_Str';
%             end
%             
%             
%             NRN(nrn_cnt).Hemisphere = Dset.SP(ii).Hemisphere;
%             NRN(nrn_cnt).APmm = Dset.SP(ii).APmm;
%             NRN(nrn_cnt).MLmm = Dset.SP(ii).MLmm;
%             if ~isempty(Dset.SP(ii).Tetrode_Channels)
%                 NRN(nrn_cnt).Tetrode_Channel_1 = Dset.SP(ii).Tetrode_Channels(1);
%             else
%                 NRN(nrn_cnt).Tetrode_Channel_1 = nan;
%             end
%             
%             NRN(nrn_cnt).Q_All = Dset.Qall(:,ii);
%             NRN(nrn_cnt).Q_Base = Dset.Qbase(:,ii);
%             NRN(nrn_cnt).Q_Post1 = Dset.Qpost1(:,ii);
%             NRN(nrn_cnt).Q_Post2 = Dset.Qpost2(:,ii);
%             NRN(nrn_cnt).Q_Post3 = Dset.Qpost3(:,ii);
%             NRN(nrn_cnt).Frate_base = Dset.FRates_base(ii);
%             NRN(nrn_cnt).Frate_post1 = Dset.FRates_post1(ii);
%             NRN(nrn_cnt).Frate_post2 = Dset.FRates_post2(ii);
%             NRN(nrn_cnt).Frate_post3 = Dset.FRates_post3(ii);
%             
%             NRN(nrn_cnt).LV_base = Dset.LocVar_base_post1_post2_post3(ii,1);
%             NRN(nrn_cnt).LV_post1 = Dset.LocVar_base_post1_post2_post3(ii,2);
%             NRN(nrn_cnt).LV_post2 = Dset.LocVar_base_post1_post2_post3(ii,3);
%             NRN(nrn_cnt).LV_post3 = Dset.LocVar_base_post1_post2_post3(ii,4);
%             
%             %         NRN(nrn_cnt).IMU_All = Dset.IMU_Speed_All;
%             %         NRN(nrn_cnt).IMU_Base = Dset.IMU_Speed_base;
%             %         NRN(nrn_cnt).IMU_Post1 = Dset.IMU_Speed_post1;
%             %         NRN(nrn_cnt).IMU_Post2 = Dset.IMU_Speed_post2;
%             %         NRN(nrn_cnt).IMU_Post3 = Dset.IMU_Speed_post3;
%             
%             NRN(nrn_cnt).Corr_speed_FR_r_base = Dset.Corr_speed_FR_r_p_base(1,ii);
%             NRN(nrn_cnt).Corr_speed_FR_r_post1 = Dset.Corr_speed_FR_r_p_post1(1,ii);
%             NRN(nrn_cnt).Corr_speed_FR_r_post2 = Dset.Corr_speed_FR_r_p_post2(1,ii);
%             NRN(nrn_cnt).Corr_speed_FR_r_post3 = Dset.Corr_speed_FR_r_p_post3(1,ii);
%             NRN(nrn_cnt).Corr_speed_FR_p_base = Dset.Corr_speed_FR_r_p_base(2,ii);
%             NRN(nrn_cnt).Corr_speed_FR_p_post1 = Dset.Corr_speed_FR_r_p_post1(2,ii);
%             NRN(nrn_cnt).Corr_speed_FR_p_post2 = Dset.Corr_speed_FR_r_p_post2(2,ii);
%             NRN(nrn_cnt).Corr_speed_FR_p_post3 = Dset.Corr_speed_FR_r_p_post3(2,ii);
%             
%             
%             NRN(nrn_cnt).AC_base = squeeze(Dset.AC_base_post1_post2_post3(ii,:,1));
%             NRN(nrn_cnt).AC_post1 = squeeze(Dset.AC_base_post1_post2_post3(ii,:,2));
%             NRN(nrn_cnt).AC_post2 = squeeze(Dset.AC_base_post1_post2_post3(ii,:,3));
%             NRN(nrn_cnt).AC_post3 = squeeze(Dset.AC_base_post1_post2_post3(ii,:,4));
%             
%             wv = Dset.SP(ii).WV_LONG.mn';
%             [~,ix] = max(max(wv));
%             ALL_WV(nrn_cnt,:) = wv(:,ix)';
%             
%             ALL_AC_base(nrn_cnt,:) = Dset.AC_base_post1_post2_post3(ii,:,1);
%             
%             nrn_cnt = nrn_cnt + 1;
%         end
%         
%      IMU_cnt = IMU_cnt + 1;
% 
%      
%     else
%     end
%     
% end

% save to table
TBL = struct2table(NRN);
% figure(150)
% for io = 1:length(ALL_WV)
%     plot(ALL_WV(io,:))
%     pause
% end

ALL_WV = ALL_WV * -1;
IX_M1 = categorical(TBL.BrainRegion) == 'M1';
M1_WV = ALL_WV(IX_M1,:);

% test = logical(M1_WV_yesno); 
% M1_WV = M1_WV(test,:);
% M1_WV = M1_WV * -1;

WV_point = (1/30000)*1000;
ALL_WV_x_msec = 0:WV_point:WV_point*57;

M1_AC_base = ALL_AC_base(IX_M1,:);
M1_AC_base = M1_AC_base(IX_good_spikes,:);
M1_WV = M1_WV(IX_good_spikes,:);
AC_x_msec = Dset.AC_x_ms;
% creating waveform features
[peak_half_width, peak_to_trough] = Spike_width(M1_WV);
wv_features = [peak_half_width, peak_to_trough];% figure; scatter(peak_half_width, peak_to_trough)
wv_feature_labels = {'half width' 'pk to tr'};
plot_it = true;

[N_type,type_lables] = neuron_subtypes(wv_features,wv_feature_labels,M1_WV,ALL_WV_x_msec,M1_AC_base,AC_x_msec,plot_it);

TBL_M1 = TBL(IX_M1,:);
TBL_M1 = TBL_M1(test,:);

TBL_M1.Neuron_type = N_type;

% Find FR rates of 0 during baseline or any drug period and exclude from analysis
IX_noFR = TBL_M1.Frate_base == 0 | TBL_M1.Frate_post1 == 0 | TBL_M1.Frate_post2 == 0 ...
           | TBL_M1.Frate_post3 == 0;
       
IX_LID = categorical(TBL_M1.Group) == '6ODHA_LID';
IX_SHAM = categorical(TBL_M1.Group) == 'SHAM';

    
corr_yes_LID = length(find(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_LID) < .05));
corr_no_LID = length(find(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_LID) > .05));
% corr_nan_LID = length(find(isnan(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_LID))));
corr_yes_SHAM = length(find(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_SHAM) < .05));
corr_no_SHAM = length(find(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_SHAM) > .05));
% corr_nan_SHAM = length(find(isnan(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_SHAM))));
figure; 
subplot(121)
pie([corr_yes_LID corr_no_LID])
title('LID Ket period mov corr')
subplot(122)
pie([corr_yes_SHAM corr_no_SHAM])
title('SHAM Ket period mov corr')
% title('Proportion of M1 cells correlated with IMU speed during ketamine period across 5 rats')

corr_yes_L_SHAM = length(find(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_SHAM) < .05));
corr_no_L_SHAM = length(find(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_SHAM) > .05));
corr_nan_L_SHAM = length(find(isnan(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_SHAM))));
corr_yes_L_LID = length(find(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_LID) < .05));
corr_no_L_LID = length(find(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_LID) > .05));
corr_nan_L_LID = length(find(isnan(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_LID))));
figure; 
subplot(121)
pie([corr_yes_L_LID corr_no_L_LID])
title('LID LDOPA period mov corr')
subplot(122)
pie([corr_yes_L_SHAM corr_no_L_SHAM])
title('SHAM LDOPA period mov corr')

% corr_yes_L = length(find(TBL_M1.Corr_speed_FR_p_post1 < .05));
% corr_no_L = length(find(TBL_M1.Corr_speed_FR_p_post1 > .05));
% corr_nan_L = length(find(isnan(TBL_M1.Corr_speed_FR_p_post1)));
% figure; pie([corr_yes_L corr_no_L corr_nan_L])
% title('Proportion of M1 cells correlated with IMU speed during peak LDOPA period across 5 rats')
% 
corr_yes_B_SHAM = length(find(TBL_M1.Corr_speed_FR_p_base(~IX_noFR & IX_SHAM) < .05));
corr_no_B_SHAM = length(find(TBL_M1.Corr_speed_FR_p_base(~IX_noFR & IX_SHAM) > .05));
corr_nan_B_SHAM = length(find(isnan(TBL_M1.Corr_speed_FR_p_base(~IX_noFR & IX_SHAM))));
corr_yes_B_LID = length(find(TBL_M1.Corr_speed_FR_p_base(~IX_noFR & IX_LID) < .05));
corr_no_B_LID = length(find(TBL_M1.Corr_speed_FR_p_base(~IX_noFR & IX_LID) > .05));
corr_nan_B_LID = length(find(isnan(TBL_M1.Corr_speed_FR_p_base(~IX_noFR & IX_LID))));
figure; 
subplot(121)
pie([corr_yes_B_LID corr_no_B_LID])
title('LID baseline period mov corr')
subplot(122)
pie([corr_yes_B_SHAM corr_no_B_SHAM])
title('SHAM baseline period mov corr')

corr_yes_LID_end = length(find(TBL_M1.Corr_speed_FR_p_post3(~IX_noFR & IX_LID) < .05));
corr_no_LID_end = length(find(TBL_M1.Corr_speed_FR_p_post3(~IX_noFR & IX_LID) > .05));
% corr_nan_LID = length(find(isnan(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_LID))));
corr_yes_SHAM_end = length(find(TBL_M1.Corr_speed_FR_p_post3(~IX_noFR & IX_SHAM) < .05));
corr_no_SHAM_end = length(find(TBL_M1.Corr_speed_FR_p_post3(~IX_noFR & IX_SHAM) > .05));

figure
bar([corr_yes_B_LID corr_no_B_LID; corr_yes_L_LID corr_no_L_LID; ...
    corr_yes_LID corr_no_LID; corr_yes_LID_end corr_no_LID_end])
pubify_figure_axis
title('LID animal neurons correlated with IMU')

figure
bar([corr_yes_B_SHAM corr_no_B_SHAM; corr_yes_L_SHAM corr_no_L_SHAM; ...
    corr_yes_SHAM corr_no_SHAM; corr_yes_SHAM_end corr_no_SHAM_end])
pubify_figure_axis
title('SHAM animal neurons correlated with IMU')

% figure; pie([corr_yes_B_SHAM corr_no_B_SHAM corr_nan_B_SHAM])
% title('Proportion of M1 cells correlated with IMU speed during baseline across 5 rats')

Corrpos_Ket_IX = sum(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_SHAM) < .05 & TBL_M1.Corr_speed_FR_r_post2(~IX_noFR & IX_SHAM) > 0);
Corrneg_Ket_IX = sum(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_SHAM) < .05 & TBL_M1.Corr_speed_FR_r_post2(~IX_noFR & IX_SHAM) < 0);
figure; pie([Corrpos_Ket_IX Corrneg_Ket_IX])
title('SHAM positively and negatively corr during ketamine')

Corrpos_Ket_IX_LID = sum(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_LID) < .05 & TBL_M1.Corr_speed_FR_r_post2(~IX_noFR & IX_LID) > 0);
Corrneg_Ket_IX_LID = sum(TBL_M1.Corr_speed_FR_p_post2(~IX_noFR & IX_LID) < .05 & TBL_M1.Corr_speed_FR_r_post2(~IX_noFR & IX_LID) < 0);
figure; pie([Corrpos_Ket_IX_LID Corrneg_Ket_IX_LID])
title('LID positively and negatively corr during ketamine')

Corrpos_LDOPA_IX = sum(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_SHAM) < .05 & TBL_M1.Corr_speed_FR_r_post1(~IX_noFR & IX_SHAM) > 0);
Corrneg_LDOPA_IX = sum(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_SHAM) < .05 & TBL_M1.Corr_speed_FR_r_post1(~IX_noFR & IX_SHAM) < 0);
figure; pie([Corrpos_LDOPA_IX Corrneg_LDOPA_IX])
title('SHAM positively and negatively corr during during peak LDOPA')

Corrpos_LDOPA_IX_LID = sum(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_LID) < .05 & TBL_M1.Corr_speed_FR_r_post1(~IX_noFR & IX_LID) > 0);
Corrneg_LDOPA_IX_LID = sum(TBL_M1.Corr_speed_FR_p_post1(~IX_noFR & IX_LID) < .05 & TBL_M1.Corr_speed_FR_r_post1(~IX_noFR & IX_LID) < 0);
figure; pie([Corrpos_LDOPA_IX_LID Corrneg_LDOPA_IX_LID])
title('LID positively and negatively corr during during peak LDOPA')

Corrpos_base_IX = sum(TBL_M1.Corr_speed_FR_p_base < .05 & TBL_M1.Corr_speed_FR_r_base > 0);
Corrneg_base_IX = sum(TBL_M1.Corr_speed_FR_p_base < .05 & TBL_M1.Corr_speed_FR_r_base < 0);
figure; pie([Corrpos_base_IX Corrneg_base_IX])
title('Sig Proportion of M1 cells positively and negatively corr with IMU speed during baseline across 5 rats')

% Firing rate changes before and after ketamine and LDOPA
TBL_M1.Frate_post2mbase = (TBL_M1.Frate_post2 - TBL_M1.Frate_base)./(TBL_M1.Frate_post2 + TBL_M1.Frate_base);
TBL_M1.Frate_post1mbase = (TBL_M1.Frate_post1 - TBL_M1.Frate_base)./(TBL_M1.Frate_post1 + TBL_M1.Frate_base);
figure; boxplot(TBL_M1.Frate_post2mbase(~IX_noFR,:), {TBL_M1.Group(~IX_noFR,:) TBL_M1.Neuron_type(~IX_noFR,:)})
title('Post ket versus baseline normalised M1 between groups and neuron type')

figure; boxplot(TBL_M1.Frate_post1mbase(~IX_noFR,:), {TBL_M1.Group(~IX_noFR,:) TBL_M1.Neuron_type(~IX_noFR,:)})
title('Post LDOPA versus baseline normalised M1 between groups and neuron type')

[h, p] = ttest(TBL_M1.Frate_post1mbase(~IX_noFR & IX_SHAM));
% Local variance difference
TBL_M1.LV_post2mbase = TBL_M1.LV_post2 - TBL_M1.LV_base;
TBL_M1.LV_post1mbase = TBL_M1.LV_post1 - TBL_M1.LV_base;

figure
histogram_cowen({TBL_M1.LV_post1mbase(~IX_noFR & IX_LID) TBL_M1.LV_post2mbase(~IX_noFR & IX_LID)}, .05)

figure
histogram_cowen({TBL_M1.LV_post1mbase(~IX_noFR & IX_SHAM) TBL_M1.LV_post2mbase(~IX_noFR & IX_SHAM)}, .05)

figure
histogram_cowen({TBL_M1.LV_base(~IX_noFR,:) TBL_M1.LV_post1(~IX_noFR,:) TBL_M1.LV_post2(~IX_noFR,:)}, .05)

