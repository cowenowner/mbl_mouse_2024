%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q11_Abhi_Does_ketamine_affect_burstiness_of_neurons';
ses_to_ana = 'Q11_Abhi_Does_ketamine_affect_burstiness_of_neurons';

GP.Analysis_Dir = 'E:\Temp\TempAnaResults';

PLOT_IT = false;
adir = fullfile(GP.Analysis_Dir,ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end
cm = lines(5);

ses = {};
ALL = [];
ALL.Control = [];
ALL.LID = [];
ALL.animal = [];
ALL.SESSIONID = [];
ALL.SES.day = [];
ALL.SES.animal = [];
ALL.SES.group = [];

%%%%%%%%%%%%%%%%%%%%
NRN = [];
ACb = [];
ACp1 = [];
ACp2 = [];

AC2b = [];
AC2p1 = [];
AC2p2 = [];

AllQbase = [];
AllQset = [];
AllQ = cell(3,1);
AllQ_all  = [];

All_TS_al = [];
All_TS_depth = [];


ALL_R = [];
ALL_R_p = [];
AllRset = [];
WV= [];
ses_cnt = 1;
iEpoch =2;
nrn_cnt = 1;

for iF = 1:length(d)

    Dset = load(fullfile(adir,d(iF).name));
    [Dset.TM.BaselineUsec/60e6 diff(Dset.TM.BaselineUsec/60e6);
        Dset.TM.EventStartEndUsec/60e6 diff(Dset.TM.EventStartEndUsec/60e6);
        Dset.TM.PostInjectionUsec(1,:)/60e6 diff(Dset.TM.PostInjectionUsec(1,:)/60e6);
        Dset.TM.PostInjectionUsec(2,:)/60e6 diff(Dset.TM.PostInjectionUsec(2,:)/60e6)
        ];
    
%     Dset.LocVar_Before_Early_Late_M1 = [];
%     Dset.LocVar_Before_Early_Late_Striatum = [];
%     Dset.FRates_base_post1_M1 = [];
%     Dset.FRates_base_post1_Striatum = [];
%     
%     for iD = 1:length(Dset.SP)
%         if Dset.SP(iD).Depth_uM < 3000 
%             Dset.LocVar_Before_Early_Late_M1(iD,:) = Dset.LocVar_Before_Early_Late(iD,:);
%             Dset.LocVar_Before_Early_Late_Striatum(iD,1:3) = nan;
%             Dset.FRates_base_post1_M1(iD,:) = [Dset.FRates_base(iD) Dset.FRates_post1(iD)];
%             Dset.FRates_base_post1_Striatum(iD,1:2) = nan;
%         else 
%             Dset.LocVar_Before_Early_Late_Striatum(iD,:) = Dset.LocVar_Before_Early_Late(iD,:);
%             Dset.LocVar_Before_Early_Late_M1(iD,1:3) = nan;
%             Dset.FRates_base_post1_Striatum(iD,:) = [Dset.FRates_base(iD) Dset.FRates_post1(iD)];
%             Dset.FRates_base_post1_M1(iD,1:2) = nan;
%         end
%     end
%     
%     for iFR = 1:length(Dset.SP)
%         if ISNAN(Dset.LocVar_Before_Early_Late(iFR,1:3))
%             Dset.FRates_base_post1_M1(iFR,:) = nan;
%             Dset.FRates_base_post1_Striatum(iFR,:) = nan;
%         else
%         end
%     end
%     
%     
%     IX_M1 = Dset.LocVar_Before_Early_Late_M1(:,1) >= 1;
%     if Dset.SES.RatType == 'Control' 
%         ALL.Control.LocVar_Before_Early_Late_M1(iF,:) = nanmean(Dset.LocVar_Before_Early_Late(IX_M1,:));
%         ALL.Control.FRates_base_post1_M1(iF,:)
%     elseif Dset.SES.RatType == '6ODHA_LID'
%         ALL.LID.LocVar_Before_Early_Late_M1(iF,:) = nanmean(Dset.LocVar_Before_Early_Late(IX_M1,:));
%     end
%     
%     if ~isempty(Dset.LocVar_Before_Early_Late_Striatum)
%         IX_Striat = Dset.LocVar_Before_Early_Late_Striatum(:,1) >= 1;
%         if Dset.SES.RatType == 'Control'
%             ALL.Control.LocVar_Before_Early_Late_Striatum(iF,:) = nanmean(Dset.LocVar_Before_Early_Late(IX_Striat,:));
%         elseif Dset.SES.RatType == '6ODHA_LID'
%             ALL.LID.LocVar_Before_Early_Late_Striatum(iF,:) = nanmean(Dset.LocVar_Before_Early_Late(IX_Striat,:));
%         end
%     end

    % the cool net level analyses.
    ALL.SES.LocVar_Before_Early_Late(iF,:) = nanmean(Dset.LocVar_Before_Early_Late);
    AllQ{1} = [AllQ{1} int16(Dset.Qbase)];
    AllQ{2} = [AllQ{2} int16(Dset.Qpost1)];
    AllQ{3} = [AllQ{3} int16(Dset.Qpost2)];
    
    AllQ_all = [AllQ_all int16(Dset.Qall)];

    AllQset = [AllQset;ones(Cols(Dset.Qpost2),1)*iF];
    % the r values.
    t = triu(ones(size(Dset.R_base)),1)==1;
    ALL_R = [ALL_R; Dset.R_base(t) Dset.R_post1(t) Dset.R_post2(t)];
    ALL_R_p = [ALL_R_p; Dset.R_p_base(t) Dset.R_p_post1(t) Dset.R_p_post2(t)];
    AllRset = [AllRset;ones(length(Dset.R_base(t)),1)*iF];
    bin_sec = (Dset.Qbase_bins_uS(2) - Dset.Qbase_bins_uS(1))/1e6;


    if PLOT_IT
        figure
        subplot(1,6,1:2)
        imagesc(mean(Dset.Qbase_bins_uS,2)/60e6,[], Dset.Qbase'/bin_sec);
        ylabel('Neuron')
        title('Baseline')
        subplot(1,6,3:4)
        imagesc(mean(Dset.Qpost1_bins_uS,2)/60e6,[], Dset.Qpost1'/bin_sec);
        title('Post 1')
        xlabel('min')
        subplot(1,6,5:6)
        imagesc(mean(Dset.Qpost2_bins_uS,2)/60e6,[], Dset.Qpost2'/bin_sec);
        title('Post 2')
        equalize_color_axes
        set(gcf,'Position',[ 41.8        503.4       1349.6        354.6])
        colorbar_label('Rate (Hz')
    end
    % Ketamine modulation
    for ii = 1:length(Dset.FRates_base)
        All_TS_al{nrn_cnt} = Dset.TSaligned{ii};

        
        NRN(nrn_cnt).NeuronID = nrn_cnt;
        if Dset.SES.RatType == 'Control'
            NRN(nrn_cnt).Group = 'Control';
        elseif Dset.SES.RatType == '6ODHA_LID'
            NRN(nrn_cnt).Group = '6ODHA_LID';
        end
        if isempty(Dset.SP(ii).Depth_uM) 
            dp = nan;
        elseif Dset.SP(ii).Depth_uM ==0
            dp = nan;
        else
            dp = Dset.SP(ii).Depth_uM;
        end
        NRN(nrn_cnt).Depth_uM = dp;
        if dp < 3000
            NRN(nrn_cnt).BrainRegion = 'M1';
        else
            NRN(nrn_cnt).BrainRegion = 'Striatum';
        end
        NRN(nrn_cnt).Hemisphere = Dset.SP(ii).Hemisphere;
        NRN(nrn_cnt).APmm = Dset.SP(ii).APmm;
        NRN(nrn_cnt).MLmm = Dset.SP(ii).MLmm;
        if ~isempty(Dset.SP(ii).Tetrode_Channels)
            NRN(nrn_cnt).Tetrode_Channel_1 = Dset.SP(ii).Tetrode_Channels(1);
        else
            NRN(nrn_cnt).Tetrode_Channel_1 = nan;
        end
        
        NRN(nrn_cnt).Session = iF;
        NRN(nrn_cnt).Frate_base = Dset.FRates_base(ii);
        NRN(nrn_cnt).Frate_post1 = Dset.FRates_post1(ii);
        NRN(nrn_cnt).Frate_post2 = Dset.FRates_post2(ii);
        NRN(nrn_cnt).Frate_post1mbase = Dset.FRates_post1(ii) -  Dset.FRates_base(ii);
        NRN(nrn_cnt).Frate_post2mbase = Dset.FRates_post2(ii) -  Dset.FRates_base(ii);
        
        NRN(nrn_cnt).LocVar_base = Dset.LocVar_Before_Early_Late(ii,1);
        NRN(nrn_cnt).LocVar_post1  = Dset.LocVar_Before_Early_Late(ii,2);
        NRN(nrn_cnt).LocVar_post2  = Dset.LocVar_Before_Early_Late(ii,3);
        NRN(nrn_cnt).LocVar_post1mbase = NRN(nrn_cnt).LocVar_post1 -  NRN(nrn_cnt).LocVar_base;
        NRN(nrn_cnt).LocVar_post2mbase = NRN(nrn_cnt).LocVar_post2 -  NRN(nrn_cnt).LocVar_base;
        
        wv = Dset.SP(ii).WV.mWV;
        [~,ix] = max(max(wv));
        WV(nrn_cnt,:) = wv(:,ix)';
        
        ACb(nrn_cnt,:) = squeeze(Dset.AC(ii,:,1));
        ACp1(nrn_cnt,:) = squeeze(Dset.AC(ii,:,2));
        ACp2(nrn_cnt,:) = squeeze(Dset.AC(ii,:,3));
        
        %         AC2b(nrn_cnt,:) = squeeze(Dset.AcorrSmth_Before_Early_Late(ii,:,1));
        %         AC2p1(nrn_cnt,:) = squeeze(Dset.AcorrSmth_Before_Early_Late(ii,:,2));
        %         AC2p2(nrn_cnt,:) = squeeze(Dset.AcorrSmth_Before_Early_Late(ii,:,3));
        %
        %         sPSDb(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_Late(ii,:,1));
        %         sPSDp1(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_Late(ii,:,2));
        %         sPSDp2(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_Late(ii,:,3));
        %
        %         sPSDbZ(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_LateZ(ii,:,1));
        %         sPSDp1Z(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_LateZ(ii,:,2));
        %         sPSDp2Z(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_LateZ(ii,:,3));
        
        nrn_cnt = nrn_cnt + 1;
    end
    
end

TBL = struct2table(NRN);
% writetable(TBL,fullfile('C:\Temp\','Q11_Abhi_ketamine_and_burst.csv'))
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the Loc Var diff
 if PLOT_IT
     figure
     subplot(311)
     histogram_cowen({ALL.LID.LocVar_Before_Early_Late_M1(:,2) - ALL.LID.LocVar_Before_Early_Late_M1(:,1)},.05)
     title('LocVar_post1mbase_M1_LID')
     subplot(312)
     histogram_cowen({ALL.LID.LocVar_Before_Early_Late_M1(:,3) - ALL.LID.LocVar_Before_Early_Late_M1(:,1)},.05)
     title('LocVar_post2mbase_M1_LID')
     subplot(313)
     histogram_cowen({ALL.LID.LocVar_Before_Early_Late_M1(:,3) - ALL.LID.LocVar_Before_Early_Late_M1(:,2)},.05)
     title('LocVar_post2post1_M1_LID')
     
     figure
     subplot(311)
     histogram_cowen({ALL.Control.LocVar_Before_Early_Late_M1(:,2) - ALL.Control.LocVar_Before_Early_Late_M1(:,1)},.05)
     title('LocVar_post1mbase_M1_Control')
     subplot(312)
     histogram_cowen({ALL.Control.LocVar_Before_Early_Late_M1(:,3) - ALL.Control.LocVar_Before_Early_Late_M1(:,1)},.05)
     title('LocVar_post2mbase_M1_Control')
     subplot(313)
     histogram_cowen({ALL.Control.LocVar_Before_Early_Late_M1(:,3) - ALL.Control.LocVar_Before_Early_Late_M1(:,2)},.05)
     title('LocVar_post2post1_M1_Control')
     
     figure
     subplot(311)
     histogram_cowen({ALL.Control.LocVar_Before_Early_Late_Striatum(:,2) - ALL.Control.LocVar_Before_Early_Late_Striatum(:,1)},.05)
     title('LocVar_post1mbase_Striatum_Control')
     subplot(312)
     histogram_cowen({ALL.Control.LocVar_Before_Early_Late_Striatum(:,3) - ALL.Control.LocVar_Before_Early_Late_Striatum(:,1)},.05)
     title('LocVar_post2mbase_Striatum_Control')
     subplot(313)
     histogram_cowen({ALL.Control.LocVar_Before_Early_Late_Striatum(:,3) - ALL.Control.LocVar_Before_Early_Late_Striatum(:,2)},.05)
     title('LocVar_post2post1_Striatum_Control')
     
 end