function OUT = Q18b_Do_spikes_phase_lock_to_80Hz_power_PETH_Abhi()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 80Hz power fluctuates a fair bit is their phase locking to the power
% envelop of 80Hz oscillations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abhi 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = [];
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();

PLOT_IT = false;
fq_ctrs = 84;
SES = LK_Session_Info();
n_ac_bins = 100;
ac_binsize_ms = 4;
ac_x_ms = ((1:n_ac_bins)*ac_binsize_ms) - ac_binsize_ms/2;

OUT.SES = SES;
OUT.aborted = false;
OUT.SP = SP;
OUT.n_ac_bins = n_ac_bins;
OUT.ac_binsize_ms = ac_binsize_ms;
tbl_cnt = 1;
DO_reref = false;
DO_regress_80 = true;
fq_ctrl_band = {[62 70]};
fq_bands = {'gamma_80'}; % [74 94] Hz
LFP_sFreq = 500;

%% Load the best non-reref.
SP_tbl = struct2table(SP);

LFPfile{1} = find_files(fullfile('Left_M1_best_nonreref_ratio*.mat'));
IX_L_M1 = SP_tbl.Tetrode < 9 & SP_tbl.Depth_uM < 2800;
IX_All(:,1) = IX_L_M1;
LFPfile{2} = find_files(fullfile('Right_M1_best_nonreref_ratio*.mat'));
IX_R_M1 = SP_tbl.Tetrode > 8 & SP_tbl.Depth_uM < 2800;
IX_All(:,2) = IX_R_M1;
cnt = 3;
if ~isempty(dir('Left_STR_best_nonreref_ratio*.mat'))
    LFPfile{cnt} = find_files(fullfile('Left_STR_best_nonreref_ratio*.mat'));
    IX_L_STR = SP_tbl.Tetrode < 9 & SP_tbl.Depth_uM > 2800;
    IX_All(:,cnt) = IX_L_STR;
    cnt = cnt +1;
else 
end

if ~isempty(dir('Right_STR_best_nonreref_ratio*.mat'))
    LFPfile{cnt} = find_files(fullfile('Right_STR_best_nonreref_ratio*.mat'));
    IX_R_STR = SP_tbl.Tetrode > 8 & SP_tbl.Depth_uM > 2800;
    IX_All(:,cnt) = IX_R_STR;
else 
end

LFP_RR = [];
if ~isempty(dir('Left_M1_700*.mat'))
    LFP_RR{1} = find_files(fullfile('Left_M1_700*.mat'));
else
    LFP_RR{1} = []; 
end
if ~isempty(dir('Right_M1_700*.mat'))
    LFP_RR{2} = find_files(fullfile('Right_M1_700*.mat'));
else
    LFP_RR{2} = [];
end


%% determine what events are available for alignment.
intervals_around_evt_min = [];
inj1_start_end_uS = [];
inj2_start_end_uS = [];

if any(E.EventID == 'LDOPAInjectionStart')
        inj1_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
        inj1_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        if any(E.EventID == 'SalineInjectionStart')
            inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'SalineInjectionStart')*60*1e6;
            inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'SalineInjectionEnd')*60*1e6;
            OUT.inj2 = 'Sal';
        elseif any(E.EventID == 'KetInjectionStart')
            inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
            inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
            OUT.inj2 = 'Ket';
        end
        intervals_around_evt_min =  [-25 -5; 2 22; 60 80; 2 22];
        big_peri_event_min = [-25 180];
        OUT.inj1 = 'Ldo';
else
    return
end

% Create some filters.
F = SPEC_create_filters(fq_bands,LFP_sFreq);
F_ctrl = SPEC_create_filters(fq_ctrl_band,LFP_sFreq);

pTBL = struct();

TIMES.Event1StartEndUsec(1,1) = inj1_start_end_uS(1);
TIMES.Event1StartEndUsec(1,2) = inj1_start_end_uS(2);
TIMES.Event2StartEndUsec(1,1) = inj2_start_end_uS(1);
TIMES.Event2StartEndUsec(1,2) = inj2_start_end_uS(2);
TIMES.IntervalUsec(1,1) = TIMES.Event1StartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.IntervalUsec(1,2) = TIMES.Event1StartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.IntervalUsec(2,1) = TIMES.Event2StartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.IntervalUsec(2,2) = TIMES.Event2StartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.IntervalUsec(3,1) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(2,1)*60*1e6;
TIMES.IntervalUsec(3,2) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(2,2)*60*1e6;
TIMES.IntervalUsec(4,1) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(3,1)*60*1e6;
TIMES.IntervalUsec(4,2) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(3,2)*60*1e6;

if length(intervals_around_evt_min(:,1)) > 3
    TIMES.IntervalUsec(5,1) = TIMES.Event1StartEndUsec(2) + intervals_around_evt_min(4,1)*60*1e6;
    TIMES.IntervalUsec(5,2) = TIMES.Event1StartEndUsec(2) + intervals_around_evt_min(4,2)*60*1e6;
else
end

TIMES.PeriEventUsec = [TIMES.Event1StartEndUsec(1) + big_peri_event_min(1)*60*1e6 TIMES.Event1StartEndUsec(1) + big_peri_event_min(2)*60*1e6 ];

OUT.TM = TIMES;
    
%%
for ix = 1:length(LFPfile)
    lfp = [];
    TSr = [];
    SPr = [];
    % Filter and upsample LFP
    filename = LFPfile{ix};
    LFP = LK_Load_and_Clean_LFP_Abhi('',filename{1});
    lfp = LFP.LFP;
    if DO_reref == true
        if ix < 3
            if ~isempty(LFP_RR{ix})
                filename2 = LFP_RR{ix};
                LFP_reref = LK_Load_and_Clean_LFP_Abhi('',filename2{1});
                OUT.LFP_rms_before_RR = rms(lfp);
                lfp = lfp - LFP_reref.LFP;
                OUT.LFP_rms_after_RR = rms(lfp);
            else
                OUT.LFP_rms_before_RR = nan;
                OUT.LFP_rms_after_RR = nan;
            end
        else
        end
    else
    end
    
    TSr = TS(IX_All(:,ix));
    SPr = SP_tbl(IX_All(:,ix),:);
    SPr = table2struct(SPr);
    for ii = 1:Rows(TIMES.IntervalUsec)
        ix1 = binsearch(LFP.t_uS , TIMES.IntervalUsec(ii,1));
        ix2 = binsearch(LFP.t_uS , TIMES.IntervalUsec(ii,2));
        L = [];
        Lp = [];
        for iFreq = 1:length(fq_bands)
            % Filter data bandpass for each freq and see if spikes phase
            % lock.
            L(:,iFreq) = filtfilt(F{iFreq},lfp(ix1:ix2));
            
            % do something crazy. Look at the response relative to power
            % in the frequency rather than the phase.
            ss = abs(hilbert(L(:,iFreq)));
            ss = ss - mean(ss);
            Lp(:,iFreq) = ss;
            
        end
                
        TSr = Restrict(TSr,LFP.t_uS(ix1),LFP.t_uS(ix2));
        
        if DO_regress_80 == true           
            iF = Cols(Lp)+1; 
            L(:,iF) = filtfilt(F_ctrl{1},lfp(ix1:ix2));
            % regress 80 from lower frequency band
            ssf = abs(hilbert(L(:,iF)));
            ss = abs(hilbert(L(:,1)));
            [~,~,ssr] = regress_cowen(ssf,ss);
            ssr = ssr - mean(ssr);
            Lp(:,iF) = ssr;
            for it = 1: length(TSr)
                [M{it}, ix, x_sec] = PETH_EEG_simple([LFP.t_uS(ix1:ix2), Lp(:,1)], TSr{14}, LFP_sFreq/2, LFP_sFreq/2,LFP_sFreq);
                [M_regress, ix, x_sec_regress] = PETH_EEG_simple([LFP.t_uS(ix1:ix2), Lp(:,1)], TSr, LFP_sFreq/2, LFP_sFreq/2,LFP_sFreq);
            end
        else
        % %
            [M, ix, x_sec] = PETH_EEG_simple([LFP.t_uS(ix1:ix2), Lp(:,1)], TSr, LFP_sFreq/2, LFP_sFreq/2,LFP_sFreq);
            M_regress = [];
            x_sec_regress = [];
            
        end

        for iN = 1:length(TSr)
            % meta stuff
            pTBL(tbl_cnt).NeuronID = SPr(iN).fname;
            %                 pTBL(tbl_cnt).UNeuronID = SES.rat + 1000*SES.session + 100000*iN;
            pTBL(tbl_cnt).Session = SES.session;
            pTBL(tbl_cnt).Rat = SES.rat;
            %                 pTBL(tbl_cnt).FreqBand = fq_bands{iFreq};
            pTBL(tbl_cnt).LFP_filename = filename{1};
            pTBL(tbl_cnt).Interval = ii;
            pTBL(tbl_cnt).Tetrode = SPr(iN).Tetrode;
            pTBL(tbl_cnt).Hemisphere = SPr(iN).Hemisphere;
            pTBL(tbl_cnt).Depth_uM = SPr(iN).Depth_uM;
            pTBL(tbl_cnt).WV = single(SPr(iN).WV_LONG.mn);
            
            %
            if length(TSr{iN})> 20
                [a,b]=AutoCorr(TSr{iN}/100,ac_binsize_ms,n_ac_bins);
            else
                a = nan(1,n_ac_bins);
                b = ac_x_ms;
            end
            pTBL(tbl_cnt).AC = a;
            pTBL(tbl_cnt).AC_x_ms = b/1000;
            %                 pTBL(tbl_cnt).LFP_uS_filt = single(L);
            pTBL(tbl_cnt).TSr = TSr(iN,1);
            
            % Calculated stuff

            pTBL(tbl_cnt).PETH_M = M;
            pTBL(tbl_cnt).PETH_x_sec = x_sec;
            pTBL(tbl_cnt).PETH_M_regress = M_regress;
            pTBL(tbl_cnt).PETH_x_sec_regress = x_sec_regress;
            
         
            tbl_cnt = tbl_cnt + 1;
        end
        
    end
end

TBL = struct2table(pTBL);

OUT.TBL = TBL;
OUT.pTBL = pTBL;
