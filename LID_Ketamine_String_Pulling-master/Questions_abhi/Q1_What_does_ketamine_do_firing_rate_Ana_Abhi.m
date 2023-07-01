%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars
mfile = 'Q1_What_does_ketamine_do_firing_rate_Abhi_Ana';
% ses_to_ana = 'Q1_What_does_ketamine_do_firing_rate_Abhi_10secbins';
% % code_dir = fileparts(which('Q1_What_does_ketamine_do_firing_rate'));
%    mfile = 'Q5_What_does_LDOPA_do_firing_rate_Ana';
%    ses_to_ana = 'Q5_What_does_LDOPA_do_firing_rate';
% code_dir = fileparts(which('Q1_What_does_ketamine_do_firing_rate'));

% GP = LK_Globals;
GP.Analysis_Dir = 'E:\Temp_6.6.23\QALL_matrix';

PLOT_IT = false;
% adir = fullfile(GP.Analysis_Dir,ses_to_ana);
adir = fullfile(GP.Analysis_Dir);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end
cm = lines(5);

ses = {};
ALL = [];
ALL.group = [];
ALL.day = [];
ALL.animal = [];
ALL.SESSIONID = [];
ALL.SES.day = [];
ALL.SES.animal = [];
ALL.SES.group = [];

%%%%%%%%%%%%%%%%%%%%
NRN = [];
ALL_AC_base = [];
ALL_WV = [];
ACb = [];
ACp1 = [];
ACp2 = [];
ACp3 = [];

AllQbase = [];
AllQset = [];
% AllQ = cell(4,1);
AllQ_all  = [];

All_TS_al = [];
All_TS_depth = [];
AllQ_b = [];
AllQ_p1 = [];
AllQ_p2 = [];
AllQ_p3 = [];

ALL_R = [];
ALL_R_p = [];
AllRset = [];
ALL_Rat = [];
ALL_IMU = [];
WV= [];
ses_cnt = 1;
iEpoch =2;
nrn_cnt = 1;
IMU_cnt = 1;
%%
for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));
    [Dset.TM.BaselineUsec/60e6 diff(Dset.TM.BaselineUsec/60e6);
%         Dset.TM.EventStartEndUsec/60e6 diff(Dset.TM.EventStartEndUsec/60e6);
        Dset.TM.PostInjectionUsec(1,:)/60e6 diff(Dset.TM.PostInjectionUsec(1,:)/60e6);
        Dset.TM.PostInjectionUsec(2,:)/60e6 diff(Dset.TM.PostInjectionUsec(2,:)/60e6);
        Dset.TM.PostInjectionUsec(3,:)/60e6 diff(Dset.TM.PostInjectionUsec(3,:)/60e6)
        ];
    if categorical(Dset.SES.Drugs) == 'Saline&Ketamine'
    
%         AllQ_all = [AllQ_all int16(Dset.Qall)];
        
%         ALL_Rat(IMU_cnt).Rat = Dset.SES.Rat;
%         if Dset.SES.RatType == 'SHAM'
%             ALL_Rat(IMU_cnt).Group = 'SHAM';
%         elseif Dset.SES.RatType == '6ODHA_LID'
%             ALL_Rat(IMU_cnt).Group = '6OHDA_LID';
%         end
%         
% %         ALL_IMU_bin(IMU_cnt,:) = Dset.IMU_All_bins_uS';
% %         ALL_IMU(IMU_cnt,:) = Dset.IMU_Speed_All';
% %         
% %         IMU_cnt = IMU_cnt + 1;
        
        ALL.SES.Speed_Base_Pre_Early_Late(iF,:) = Dset.Speed_Base_Pre_Early_Late;
        ALL.SES.Jerk_base_pre_early_late(iF,:) = Dset.Jerk_base_pre_early_late;
        ALL.SES.JerkPC1_base_pre_early_late(iF,:) = Dset.JerkPC1_base_pre_early_late;
        speed_change = Dset.Speed_Base_Pre_Early_Late(2)-Dset.Speed_Base_Pre_Early_Late(1);
        jerk_change = Dset.Jerk_base_pre_early_late(2)-Dset.Jerk_base_pre_early_late(1);
%         
% %         
% % %         the cool net level analyses.
% %             AllQ = [];
        ALL.SES.LocVar_base_pre_early_late(iF,:) = nanmean(Dset.LocVar_base_pre_early_late);
%             AllQ{1} = [AllQ{1} int16(Dset.Qbase)];
%             AllQ{2} = [AllQ{2} int16(Dset.Qpost1)];
%             AllQ{3} = [AllQ{3} int16(Dset.Qpost2)];
%             AllQ{4} = [AllQ{4} int16(Dset.Qpost3)];
            AllQ_b = [AllQ_b int16(Dset.Qbase)];
            AllQ_p1 = [AllQ_p1 int16(Dset.Qpost1)];
            AllQ_p2 = [AllQ_p2 int16(Dset.Qpost2)];
            AllQ_p3 = [AllQ_p3 int16(Dset.Qpost3)];
% %         
% %            
% %         
% %             AllQset = [AllQset;ones(Cols(Dset.Qpost2),1)*iF];
% % %         the r values.
        t = triu(ones(size(Dset.R_base)),1)==1;
        ALL_R = [ALL_R; Dset.R_base(t) Dset.R_post1(t) Dset.R_post2(t) Dset.R_post3(t)];
        ALL_R_p = [ALL_R_p; Dset.R_p_base(t) Dset.R_p_post1(t) Dset.R_p_post2(t) Dset.R_p_post3(t)];
        AllRset = [AllRset;ones(length(Dset.R_base(t)),1)*iF];
        bin_sec = (Dset.Qbase_bins_uS(2) - Dset.Qbase_bins_uS(1))/1e6;
        
%         
%         if PLOT_IT
%             Plot some session stuff
%             figure
%             subplot(1,6,1:2)
%             imagesc(mean(Dset.Qbase_bins_uS,2)/60e6,[], Dset.Qbase'/bin_sec);
%             ylabel('Neuron')
%             title('Baseline')
%             subplot(1,6,3:4)
%             imagesc(mean(Dset.Qpost1_bins_uS,2)/60e6,[], Dset.Qpost1'/bin_sec);
%             title('Post 1')
%             xlabel('min')
%             subplot(1,6,5:6)
%             imagesc(mean(Dset.Qpost2_bins_uS,2)/60e6,[], Dset.Qpost2'/bin_sec);
%             title('Post 2')
%             equalize_color_axes
%             set(gcf,'Position',[ 41.8        503.4       1349.6        354.6])
%             colorbar_label('Rate (Hz')
%         end
% %         Ketamine modulation
        for ii = 1:length(Dset.FRates_base)
            All_TS_al{nrn_cnt} = Dset.TSaligned{ii};
            
            Inj_uS = Dset.TM.Event2StartEndUsec(1,2);
            %             if categorical(Dset.SES.Drugs) == 'Saline&Ketamine'
            %                 NRN(nrn_cnt).Inj_uS = Inj_uS;
            %             else
            Inj_uS = Inj_uS+(61*60e6);
            NRN(nrn_cnt).Inj_uS = Inj_uS;
            %             end
            
            
            NRN(nrn_cnt).NeuronID = Dset.SP(ii).fname;
            NRN(nrn_cnt).Rat = Dset.SES.Rat;
            NRN(nrn_cnt).Session = Dset.SES.Session;
            if Dset.SES.RatType == 'SHAM'
                NRN(nrn_cnt).Group = 'SHAM';
            elseif Dset.SES.RatType == '6ODHA_LID'
                NRN(nrn_cnt).Group = '6OHDA_LID';
            end
            NRN(nrn_cnt).Drugs = categorical(Dset.SES.Drugs);
            if isempty(Dset.SP(ii).Depth_uM)
                dp = nan;
            elseif Dset.SP(ii).Depth_uM ==0
                dp = nan;
            else
                dp = Dset.SP(ii).Depth_uM;
            end
            
            NRN(nrn_cnt).Depth_uM = dp;
            if dp < 2800
                NRN(nrn_cnt).BrainRegion = 'M1';
            else
                NRN(nrn_cnt).BrainRegion = 'Striatum';
                NRN(nrn_cnt).M1Layer = 'Str';
            end
            
            if dp < 500
                NRN(nrn_cnt).M1Layer = '123';
            elseif dp >= 500 && dp < 900
                NRN(nrn_cnt).M1Layer = '4';
            elseif dp >= 900 && dp < 2000
                NRN(nrn_cnt).M1Layer = '56';
            elseif dp >= 2000 && dp < 2800
                NRN(nrn_cnt).M1Layer = 'M1_Str';
            end
            
            NRN(nrn_cnt).Hemisphere = Dset.SP(ii).Hemisphere;
            NRN(nrn_cnt).APmm = Dset.SP(ii).APmm;
            NRN(nrn_cnt).MLmm = Dset.SP(ii).MLmm;
            if ~isempty(Dset.SP(ii).Tetrode_Channels)
                NRN(nrn_cnt).Tetrode_Channel_1 = Dset.SP(ii).Tetrode_Channels(1);
            else
                NRN(nrn_cnt).Tetrode_Channel_1 = nan;
            end
            
            
            NRN(nrn_cnt).Speed_change = speed_change;
            NRN(nrn_cnt).Jerk_change = jerk_change;
            NRN(nrn_cnt).Q_All = Dset.Qall(:,ii);
            NRN(nrn_cnt).Q_Base = Dset.Qbase(:,ii);
            NRN(nrn_cnt).Q_Post1 = Dset.Qpost1(:,ii);
            NRN(nrn_cnt).Q_Post2 = Dset.Qpost2(:,ii);
            NRN(nrn_cnt).Q_Post3 = Dset.Qpost3(:,ii);
            NRN(nrn_cnt).Frate_base = Dset.FRates_base(ii);
            NRN(nrn_cnt).Frate_post1 = Dset.FRates_post1(ii);
            NRN(nrn_cnt).Frate_post2 = Dset.FRates_post2(ii);
            NRN(nrn_cnt).Frate_post3 = Dset.FRates_post3(ii);
            NRN(nrn_cnt).Frate_post1mbase = Dset.FRates_post1(ii) -  Dset.FRates_base(ii);
            NRN(nrn_cnt).Frate_post2mbase = Dset.FRates_post2(ii) -  Dset.FRates_base(ii);
            NRN(nrn_cnt).Frate_post2mpost1 = Dset.FRates_post2(ii) -  Dset.FRates_post1(ii);
            NRN(nrn_cnt).Frate_post3mbase = Dset.FRates_post3(ii) -  Dset.FRates_base(ii);
            NRN(nrn_cnt).Frate_post3mpost1 = Dset.FRates_post3(ii) -  Dset.FRates_post1(ii);
            NRN(nrn_cnt).Frate_post1mbasezscore = (Dset.FRates_post1(ii) -  Dset.FRates_base(ii))/(Dset.FRates_post1(ii) + Dset.FRates_base(ii));
            NRN(nrn_cnt).Frate_post2mbasezscore = (Dset.FRates_post2(ii) -  Dset.FRates_base(ii))/(Dset.FRates_post2(ii) + Dset.FRates_base(ii));
            NRN(nrn_cnt).Frate_post2mpost1zscore = (Dset.FRates_post2(ii) -  Dset.FRates_post1(ii))/(Dset.FRates_post2(ii) + Dset.FRates_post1(ii));
            NRN(nrn_cnt).Frate_post3mpost1zscore = (Dset.FRates_post3(ii) -  Dset.FRates_post1(ii))/(Dset.FRates_post3(ii) + Dset.FRates_post1(ii));
            
            NRN(nrn_cnt).LocVar_base = Dset.LocVar_base_pre_early_late(ii,1);
            NRN(nrn_cnt).LocVar_post1  = Dset.LocVar_base_pre_early_late(ii,2);
            NRN(nrn_cnt).LocVar_post2  = Dset.LocVar_base_pre_early_late(ii,3);
            NRN(nrn_cnt).LocVar_post3  = Dset.LocVar_base_pre_early_late(ii,4);
            NRN(nrn_cnt).LocVar_post1mbase = NRN(nrn_cnt).LocVar_post1 -  NRN(nrn_cnt).LocVar_base;
            NRN(nrn_cnt).LocVar_post2mbase = NRN(nrn_cnt).LocVar_post2 -  NRN(nrn_cnt).LocVar_base;
            NRN(nrn_cnt).LocVar_post2mpost1 = NRN(nrn_cnt).LocVar_post2 -  NRN(nrn_cnt).LocVar_post1;
            NRN(nrn_cnt).LocVar_post3mbase = NRN(nrn_cnt).LocVar_post3 -  NRN(nrn_cnt).LocVar_base;
            NRN(nrn_cnt).LocVar_post3mpost1 = NRN(nrn_cnt).LocVar_post3 -  NRN(nrn_cnt).LocVar_post1;
            NRN(nrn_cnt).LocVar_post1mbasezscore = (NRN(nrn_cnt).LocVar_post1 -  NRN(nrn_cnt).LocVar_base) / (NRN(nrn_cnt).LocVar_post1 +  NRN(nrn_cnt).LocVar_base);
            NRN(nrn_cnt).LocVar_post2mbasezscore = (NRN(nrn_cnt).LocVar_post2 -  NRN(nrn_cnt).LocVar_base) / (NRN(nrn_cnt).LocVar_post2 +  NRN(nrn_cnt).LocVar_base);
            NRN(nrn_cnt).LocVar_post2mpost1zscore = (NRN(nrn_cnt).LocVar_post2 -  NRN(nrn_cnt).LocVar_post1) / (NRN(nrn_cnt).LocVar_post2 +  NRN(nrn_cnt).LocVar_post1);
            NRN(nrn_cnt).LocVar_post3mpost1zscore = (NRN(nrn_cnt).LocVar_post3 -  NRN(nrn_cnt).LocVar_post1) / (NRN(nrn_cnt).LocVar_post3 +  NRN(nrn_cnt).LocVar_post1);
            
            NRN(nrn_cnt).AC_base = squeeze(Dset.AC_base_pre_early_late(ii,:,1));
            NRN(nrn_cnt).AC_post1 = squeeze(Dset.AC_base_pre_early_late(ii,:,2));
            NRN(nrn_cnt).AC_post2 = squeeze(Dset.AC_base_pre_early_late(ii,:,3));
            NRN(nrn_cnt).AC_post3 = squeeze(Dset.AC_base_pre_early_late(ii,:,4));
            
            wv = Dset.SP(ii).WV_LONG.mn';
            [~,ix] = max(max(wv));
            ALL_WV(nrn_cnt,:) = wv(:,ix)';
            
            ALL_AC_base(nrn_cnt,:) = Dset.AC_base_pre_early_late(ii,:,1);
            ACb(nrn_cnt,:) = Dset.AC_base_pre_early_late(ii,:,1);
            ACp1(nrn_cnt,:) = Dset.AC_base_pre_early_late(ii,:,2);
            ACp2(nrn_cnt,:) = Dset.AC_base_pre_early_late(ii,:,3);
            ACp3(nrn_cnt,:) = Dset.AC_base_pre_early_late(ii,:,4);
           
            
            nrn_cnt = nrn_cnt + 1;
        end
    else
    end
end

%% Save it to a nice table that is easy to read in R.
TBL = struct2table(NRN);

% classify cell types
% define times for AC and WV
AC_x_msec = Dset.AC_x_ms; 

ALL_WV = ALL_WV * -1;

WV_point = (1/30000)*1000;
ALL_WV_x_msec = 0:WV_point:WV_point*57;
% creating waveform features
[peak_half_width, peak_to_trough] = Spike_width(ALL_WV);
wv_features = [peak_half_width, peak_to_trough];%
wv_feature_labels = {'half width' 'pk to tr'};
% classify subtypes in M1
IX_M1 = categorical(TBL.BrainRegion) == 'M1';
IX_Str = categorical(TBL.BrainRegion) == 'Striatum';

plot_it = true;
% ALL_AC_base = TBL.AC_base(:,:);
[M1_type,type_lables] = neuron_subtypes(wv_features(IX_M1,:),wv_feature_labels,ALL_WV(IX_M1,:),ALL_WV_x_msec,ALL_AC_base(IX_M1,:),AC_x_msec,plot_it, [7 11]);
[Str_type,type_lables] = neuron_subtypes(wv_features(IX_Str,:),wv_feature_labels,ALL_WV(IX_Str,:),ALL_WV_x_msec,ALL_AC_base(IX_Str,:),AC_x_msec,plot_it, [6 10]);

TBL_M1 = TBL(IX_M1,:);
TBL_Str = TBL(IX_Str,:);
TBL_M1.NeuronType = M1_type;
TBL_Str.NeuronType = Str_type;

writetable(TBL_M1,'TBL_M1.csv')
writetable(TBL_Str,'TBL_Str.csv')

for ii = 1:length(TBL.Hemisphere)
    if strcmp(TBL.Hemisphere(ii),'R') || strcmp(TBL.Hemisphere(ii),'L')
    else
        TBL.Hemisphere(ii) = 'U';
    end
end

for ii = 1:size(TBL)
   idx = find(strcmp(TBL.NeuronID{ii}, TBL_M1.NeuronID));
   if ~isempty(idx)
       TBL.Neuron_type(ii) = TBL_M1.NeuronType(idx);
   end
end
for ii = 1:size(TBL)
   idx = find(strcmp(TBL.NeuronID{ii}, TBL_Str.NeuronID));
   if ~isempty(idx)
       TBL.Neuron_type(ii) = TBL_Str.NeuronType(idx);
   end
end
%% For plotting in R
vars = {'NeuronID','Rat','Group','Drugs','Hemisphere','NeuronType','Frate_post2mpost1zscore','LocVar_post2mpost1zscore'};
M1_FR_LV = TBL_M1(:,vars);
writetable(M1_FR_LV,'M1_FR_LV.csv')
% TBL = struct2table(NRN);
% writetable(TBL,fullfile('C:\Temp\', [mfile '.csv']))
% Save all results to a monster file for easier post analysis.
% LBLS = {'base','post 1', 'post 2'};
% save(fullfile('C:\Temp\',mfile))
% load(fullfile('C:\Temp\',mfile))

%% Behavior
% makes animals move.
Q1_Behavior_Ana_sub

%% Entire population analysis.
% Firing rates go up, but some go down.
Q1_Pop_and_Rate_Ana_sub
Q1_Pop_and_Rate_Ana_sub2
%% Autocorr analysis... Also includes local variance
% Upshot: LV reduces a lot and the acorr reflects this. Kills bursting.
Q1_ACorr_and_LV_Ana_sub
%% XCorr analysis... Also includes local variance


