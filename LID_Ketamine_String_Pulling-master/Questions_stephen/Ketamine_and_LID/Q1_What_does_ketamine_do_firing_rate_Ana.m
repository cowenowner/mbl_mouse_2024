%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q1_What_does_ketamine_do_firing_rate_Ana';
ses_to_ana = 'Q1_What_does_ketamine_do_firing_rate';
% % code_dir = fileparts(which('Q1_What_does_ketamine_do_firing_rate'));
%    mfile = 'Q5_What_does_LDOPA_do_firing_rate_Ana';
%    ses_to_ana = 'Q5_What_does_LDOPA_do_firing_rate';
% code_dir = fileparts(which('Q1_What_does_ketamine_do_firing_rate'));

GP = LK_Globals;
GP.Analysis_Dir = Analysis_dir;

PLOT_IT = false;
adir = fullfile(GP.Analysis_Dir,ses_to_ana);
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
ACb = [];
ACp1 = [];
ACp2 = [];

AC2b = [];
AC2p1 = [];
AC2p2 = [];


sPSDb = [];
sPSDp1 = [];
sPSDp2 = [];

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
        ]
    ALL.SES.Speed_Before_Early_Late(iF,:) = Dset.Speed_Before_Early_Late;
    ALL.SES.Jerk_Before_Early_Late(iF,:) = Dset.Jerk_Before_Early_Late;
    ALL.SES.JerkPC1_Before_Early_Late(iF,:) = Dset.JerkPC1_Before_Early_Late;
    speed_change = Dset.Speed_Before_Early_Late(2)-Dset.Speed_Before_Early_Late(1);
    jerk_change = Dset.Jerk_Before_Early_Late(2)-Dset.Jerk_Before_Early_Late(1);

    
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
        % Plot some session stuff
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
        if isempty(Dset.SP(ii).Depth_uM) 
            dp = nan;
        elseif Dset.SP(ii).Depth_uM ==0
            dp = nan;
        else
            dp = Dset.SP(ii).Depth_uM;
        end
        NRN(nrn_cnt).Depth_uM = dp;
        NRN(nrn_cnt).Hemisphere = Dset.SP(ii).Hemisphere;
        NRN(nrn_cnt).APmm = Dset.SP(ii).APmm;
        NRN(nrn_cnt).MLmm = Dset.SP(ii).MLmm;
        if ~isempty(Dset.SP(ii).Tetrode_Channels)
            NRN(nrn_cnt).Tetrode_Channel_1 = Dset.SP(ii).Tetrode_Channels(1);
        else
            NRN(nrn_cnt).Tetrode_Channel_1 = nan;
        end
        
        NRN(nrn_cnt).Session = iF;
        NRN(nrn_cnt).Speed_change = speed_change;
        NRN(nrn_cnt).Jerk_change = jerk_change;
        NRN(nrn_cnt).Frate_base = Dset.FRates_base(ii);
        NRN(nrn_cnt).Frate_post1 = Dset.FRates_post1(ii);
        NRN(nrn_cnt).Frate_post2 = Dset.FRates_post2(ii);
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

% Save it to a nice table that is easy to read in R.
TBL = struct2table(NRN);
for ii = 1:length(TBL.Hemisphere)
    if strcmp(TBL.Hemisphere(ii),'R') || strcmp(TBL.Hemisphere(ii),'L')
    else
        TBL.Hemisphere(ii) = 'U';
    end
end
try
    % Frustrating, - sometimes works, sometimes not. No clue why.
    TBL.Hemisphere = categorical(TBL.Hemisphere);
end
writetable(TBL,fullfile('C:\Temp\', [mfile '.csv']))
% Save all results to a monster file for easier post analysis.
LBLS = {'base','post 1', 'post 2'};
save(fullfile('C:\Temp\',mfile))
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


