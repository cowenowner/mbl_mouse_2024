function OUT = Q19_Is_there_cross_freq_coupling_LID_ket_Abhi()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFC during LID 80 Hz and ketamine gamma 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abhi 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = [];
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();


DO_reref = false;
SES = LK_Session_Info();
OUT.SES = SES;
OUT.aborted = false;
OUT.SP = SP;
PLOT_IT = false;
tbl_cnt = 1;
fqs = 1:.5:170;
LFP_sFreq = 500;
low_range = [2:.2:12];
method = 'tort'; % using Tort to compare with Petterson group
DO_regress_80 = true;
fq_ctrl_band = {[62 70]};

fq_bands = {[1 4] 'theta' 'gamma_80'}; 

%% Load the best non-reref.

LFPfile{1} = find_files(fullfile('Left_M1_best_nonreref_ratio*.mat'));

LFPfile{2} = find_files(fullfile('Right_M1_best_nonreref_ratio*.mat'));
cnt = 3;
if ~isempty(dir('Left_STR_best_nonreref_ratio*.mat'))
    LFPfile{cnt} = find_files(fullfile('Left_STR_best_nonreref_ratio*.mat'));
    cnt = cnt + 1;
else 
end

if ~isempty(dir('Right_STR_best_nonreref_ratio*.mat'))
    LFPfile{cnt} = find_files(fullfile('Right_STR_best_nonreref_ratio*.mat'));
else 
end

ch_num = Numbers_from_filename_Abhi(LFPfile);

% Load LFP to tetrode conversion excel
LtT = readtable('../LFP_to_Tetrode.xlsx');
TT_num = [];
for ii = 1:length(ch_num)
    IX = LtT.LFPChannel == ch_num(ii,2);
    if any(IX)
        TT_num(ii) = LtT.Tetrode(IX);
    end
end
% get depths of lfp TT
Depth_lf = [];
for ii = 1:length(TT_num)
    IX = DEPTHS(:,1) == TT_num(ii);
    if any(IX)
        Depth_lf(ii) = DEPTHS(IX,2);
    end
end

% Load reref LFP filenames
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


if any(E.EventID == 'KetInjectionStart') && any(E.EventID == 'SalineInjectionStart')
        inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
        inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
        inj1_start_end_uS(1) = E.MinFromStart(E.EventID == 'SalineInjectionStart')*60*1e6;
        inj1_start_end_uS(2) = E.MinFromStart(E.EventID == 'SalineInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 2 22; 60 80];
        big_peri_event_min = [-25 110];
        OUT.inj1 = 'Sal';
        OUT.inj2 = 'Ket';

elseif any(E.EventID == 'LDOPAInjectionStart')
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
end

% Create some filters.
F = SPEC_create_filters(fq_bands,LFP_sFreq);
F_ctrl = SPEC_create_filters(fq_ctrl_band,LFP_sFreq);

cTBL = struct();

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

    for ii = 1:Rows(TIMES.IntervalUsec)
        ix1 = binsearch(LFP.t_uS , TIMES.IntervalUsec(ii,1));
        ix2 = binsearch(LFP.t_uS , TIMES.IntervalUsec(ii,2));
        signal = [];
        signal = lfp(ix1:ix2);
        
        CFC = SPEC_cross_fq_coupling_comod_dupre2017_Abhi(signal, LFP_sFreq, low_range, method);
        
        Lp = [];
        for iFreq = 1:length(fq_bands)
            % Filter data bandpass for each freq and see if spikes phase
            % lock.
            L(:,iFreq) = filtfilt(F{iFreq},lfp(ix1:ix2));
            
            % do something crazy. Look at the response relative to power
            % in the frequency rather than the phase.
            ss = abs(hilbert(L(:,iFreq)));
            %ss = ss - mean(ss);
            Lp(:,iFreq) = ss;
            
        end
        t_uS = LFP.t_uS(ix1:ix2);
        
        [delta_pk_idx,~,~,~] = Find_peaks_troughs_zeros(L(:,1)); 
        [theta_pk_idx,~,~,~] = Find_peaks_troughs_zeros(L(:,2)); 
        
         if DO_regress_80 == true           
            iF = Cols(Lp)+1; 
            L(:,iF) = filtfilt(F_ctrl{1},lfp(ix1:ix2));
            % regress 80 from lower frequency band
            ssf = abs(hilbert(L(:,iF)));
            ss = abs(hilbert(L(:,3)));
            [~,~,ssr] = regress_cowen(ssf,ss);
            %ssr = ssr - mean(ssr);
            Lp(:,iF) = ssr;
            
            [M_delta_regress, ~, x_sec_delta_regress] = PETH_EEG_simple([t_uS Lp(:,iF)], t_uS(delta_pk_idx), LFP_sFreq/2,LFP_sFreq/2,LFP_sFreq,PLOT_IT);
            [M_theta_regress, ~, x_sec_theta_regress] = PETH_EEG_simple([t_uS Lp(:,iF)], t_uS(theta_pk_idx), LFP_sFreq/2,LFP_sFreq/2,LFP_sFreq,PLOT_IT);
        
            
         else
            
         end      

        [M_delta, ~, x_sec_delta] = PETH_EEG_simple([t_uS Lp(:,3)], t_uS(delta_pk_idx), LFP_sFreq/2,LFP_sFreq/2,LFP_sFreq,PLOT_IT);
        [M_theta, ~, x_sec_theta] = PETH_EEG_simple([t_uS Lp(:,3)], t_uS(theta_pk_idx), LFP_sFreq/2,LFP_sFreq/2,LFP_sFreq,PLOT_IT);
        
        pxx = pwelch(signal,LFP_sFreq,LFP_sFreq/2,fqs,LFP_sFreq);
        
        cTBL(tbl_cnt).Rat = SES.rat;
        cTBL(tbl_cnt).Session = SES.session;
        cTBL(tbl_cnt).Interval = ii;
        cTBL(tbl_cnt).LFP_filename = filename{1}; 
        cTBL(tbl_cnt).Depth = Depth_lf(ix);
        cTBL(tbl_cnt).TT = TT_num(ix);
        
        cTBL(tbl_cnt).CM = CFC.CM;
        cTBL(tbl_cnt).low_fq_range = CFC.low_fq_range;
        cTBL(tbl_cnt).high_fq_range = CFC.high_fq_range;
        cTBL(tbl_cnt).method = CFC.method;
        
        cTBL(tbl_cnt).fqs_psd = fqs;
        cTBL(tbl_cnt).PSD = 10*log10(pxx);
        
        cTBL(tbl_cnt).M_delta_peaks = M_delta;
        cTBL(tbl_cnt).xsec_delta_peaks = xsec_delta;
        
        cTBL(tbl_cnt).M_delta_pk_reg = M_delta_regress;
        cTBL(tbl_cnt).xsec_delta_pk_reg = xsec_delta_regress;
        
        cTBL(tbl_cnt).M_theta_peaks = M_theta;
        cTBL(tbl_cnt).xsec_theta_peaks = xsec_theta;
        
        cTBL(tbl_cnt).M_theta_pk_reg = M_theta_regress;
        cTBL(tbl_cnt).xsec_theta_pk_reg = xsec_theta_regress;
        
        
        tbl_cnt = tbl_cnt + 1;
        
    end
end

TBL = struct2table(cTBL);

OUT.TBL = TBL;
OUT.cTBL = cTBL;