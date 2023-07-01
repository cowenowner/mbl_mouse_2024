function OUT = Q16_Autocorr_CrossCorr_LFP_before_after_Ketamine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abhi 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OUT = [];

PLOT_IT = false;
LFP_epoch_s = 30;
LFP_corr_lag_s = 2;
ket_gamma_filt = [35 75];

[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;

OUT.SP = SP;
OUT.DEPTHS = DEPTHS; % Send this info out of the function for meta analysis.
OUT.META = META;
OUT.EVT = EVT;

OUT.LFP_epoch_s = LFP_epoch_s;
OUT.LFP_corr_lag_s = LFP_corr_lag_s;

fq_bands = {'delta' 'theta' 'beta' 'gamma_50' 'gamma_80'};
LFP_sFreq = 500;

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

aTBL = struct();

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


% Create some filters.
F = SPEC_create_filters(fq_bands,LFP_sFreq);

for ix = 1:length(LFPfile)

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
        L = [];
        for iFreq = 1:length(fq_bands)
            % Filter data bandpass for each freq and see if spikes phase
            % lock.
            L(:,iFreq) = filtfilt(F{iFreq},lfp(ix1:ix2));           
        end
        
        %%
        Corr_epochs = L(1,1):LFP_epoch_s*1e6:L(end,1);
        LFP_epochs = [];
        epoch_cnt = 1;
        for ii = 1:length(Corr_epochs)-1
            IX_epochs = L(:,1) >= Corr_epochs(ii) & L(:,1) < Corr_epochs(ii+1);
            if sum(IX_epochs) > 40 % this is needed because around injection times there is no LFP
                LFP_epochs(:,epoch_cnt) = L(IX_epochs,2);
                epoch_cnt = epoch_cnt + 1;
            else
            end
        end
        
        [Ket_xcorr_epochs, lags] = xcorr(LFP_epochs,LFP_corr_lag_s*LFP.sFreq,'coeff');
        Ket_acorr = [];
        multiple = 111;
        for ix = 1:111
            if ix == 1
                Ket_acorr(:,ix) = Ket_xcorr_epochs(:,ix);
            else
                Ket_acorr(:,ix) = Ket_xcorr_epochs(:,ix+(multiple*(ix-1)));
            end
        end
        time_sec_1 = -1770:30:-150;
        time_sec_2 = 150:30:1800;
        time_sec = [time_sec_1 time_sec_2];
        IX_xcorr = single(Ket_xcorr_epochs(3,:)) ~= 1;
        
        figure
        imagesc([-2 -1 0 1 2], time_sec, Ket_xcorr_epochs(:,IX_xcorr)')
        
        %%
        aTBL(tbl_cnt).Rat = SES.rat;
        aTBL(tbl_cnt).Session = SES.session;
        aTBL(tbl_cnt).Interval = ii;
        aTBL(tbl_cnt).LFP_filename = filename{1}; 
        aTBL(tbl_cnt).Depth = Depth_lf(ix);
        aTBL(tbl_cnt).TT = TT_num(ix);
        
        aTBL(tbl_cnt).CM = CFC.CM;
        aTBL(tbl_cnt).low_fq_range = CFC.low_fq_range;
        aTBL(tbl_cnt).high_fq_range = CFC.high_fq_range;
        aTBL(tbl_cnt).method = CFC.method;
        
        tbl_cnt = tbl_cnt + 1;
        
    end
end