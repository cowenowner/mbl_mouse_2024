function OUT = Q18_Do_spikes_phase_lock_to_80Hz_power_Abhi()
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
if length(SP) > 1
    SP_tbl = struct2table(SP);
else
    OUT.aborted = true;
    return
end

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
    OUT.aborted = true;
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
    
    TSR = TS(IX_All(:,ix));
    SPr = SP_tbl(IX_All(:,ix),:);
    SPr = table2struct(SPr);
    for ii = 1:Rows(TIMES.IntervalUsec)
        ix1 = binsearch(LFP.t_uS , TIMES.IntervalUsec(ii,1));
        ix2 = binsearch(LFP.t_uS , TIMES.IntervalUsec(ii,2));
        TSr = [];
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
                
        TSr = Restrict(TSR,LFP.t_uS(ix1),LFP.t_uS(ix2));
        
        if DO_regress_80 == true           
            iF = Cols(Lp)+1; 
            L(:,iF) = filtfilt(F_ctrl{1},lfp(ix1:ix2));
            % regress 80 from lower frequency band
            ssf = abs(hilbert(L(:,iF)));
            ss = abs(hilbert(L(:,1)));
            [~,~,ssr] = regress_cowen(ssf,ss);
            ssr = ssr - mean(ssr);
            Lp(:,iF) = ssr;
            
            [PowPhCr]  = SPEC_spike_field_coupling_Abhi(TSr,[LFP.t_uS(ix1:ix2), Lp],[fq_ctrs 80]);
            
        else
        % %
            [PowPhCr] = SPEC_spike_field_coupling_Abhi(TSr,[LFP.t_uS(ix1:ix2), Lp],fq_ctrs);
            
        end
        %
        %[PhC,phinfo] = SPEC_spike_field_coupling(TSr,L,fq_ctrs,'thresh_prctile_of_power',[80 10]);
        %       [PhC,phinfo] = SPEC_spike_field_coupling(TS,[LFP.t_uS(ix1:ix2), L],[8 22 50 84],'thresh_prctile_of_power',80);
        
        
        
        %         for iu = 1:length(TSr)
        %             if PhC(iu).Ang_p(4) < .05
        %                 L_res = Restrict(L,PhC(iu).good_pow_intervals{4,1});
        %                 figure;
        %                 [PhaseLocking_80,ix,x_sec] = PETH_EEG_simple([L_res(:,1), L_res(:,5)], TSr{1,iu},LFP.sFreq/50,LFP.sFreq/50,LFP.sFreq,true);
        %                 colorbar
        %                 title(sprintf ('Neuron %0.1f Rat 342 Ses 6 Peak 80Hz period', iu))
        %             end
        %         end
        
        
        %        [PhCr]  = SPEC_spike_field_coupling(TS,[LFP.t_uS(ix1:ix2), L]);
        %
        %
        %         % spike field coherence cowen
        %         [PCoh,PCohSh] = SPEC_spike_field_coherence_cowen(TS, [LFP.t_uS(ix1:ix2) LFP.LFP(ix1:ix2)], LFP.sFreq, psd_fqs);
        %
        % spike filed coupling for peak frequency
        %         [PhPk,PhPkSh] = SPEC_spike_field_coherence_track_pk_fq(TS, [LFP.t_uS(ix1:ix2) LFP.LFP(ix1:ix2)], LFP.sFreq, fqs_peakcoupling);
        
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
            pTBL(tbl_cnt).TSr = TSr{iN};
            
            % Calculated stuff
            % %                 pTBL(tbl_cnt).Ang = PhC(iN).Ang;
%             pTBL(tbl_cnt).fq_ctrs = PhC(iN).fq_ctrs;
%             pTBL(tbl_cnt).hist_rad = PhC(iN).hist_rad;
%             % pTBL(tbl_cnt).Ang_ts = PhC(iN).Ang_ts;
%             pTBL(tbl_cnt).Ang_p = PhC(iN).Ang_p;
%             pTBL(tbl_cnt).Ang_to_shuf_p = PhC(iN).Ang_to_shuf_p;
%             pTBL(tbl_cnt).Ang_to_shuf_F = PhC(iN).Ang_to_shuf_F;
%             pTBL(tbl_cnt).Ang_z = PhC(iN).Ang_z;
%             pTBL(tbl_cnt).sh_hist_rad_mn = PhC(iN).sh_hist_rad_mn;
%             pTBL(tbl_cnt).sh_hist_rad_95ci = PhC(iN).sh_hist_rad_95ci;
%             pTBL(tbl_cnt).sig_ph_locking = PhC(iN).sig_ph_locking;
%             pTBL(tbl_cnt).n_bins_above_shuff = PhC(iN).n_bins_above_shuff;
%             pTBL(tbl_cnt).good_pow_intervals = PhC(iN).good_pow_intervals;
%             
            
            %
            % %                 pTBL(tbl_cnt).rAng = PhCr(iN).Ang;
            % %                 pTBL(tbl_cnt).rhist_rad = PhCr(iN).hist_rad;
            % %                 pTBL(tbl_cnt).rAng_p = PhCr(iN).Ang_p;
            % %                 pTBL(tbl_cnt).rAng_to_shuf_p = PhCr(iN).Ang_to_shuf_p;
            % %                 pTBL(tbl_cnt).rAng_z = PhCr(iN).Ang_z;
            % %                 pTBL(tbl_cnt).rsh_hist_rad_mn = PhCr(iN).sh_hist_rad_mn;
            % %                 pTBL(tbl_cnt).rsh_hist_rad_95ci = PhCr(iN).sh_hist_rad_95ci;
            % %
            %
            %                 pTBL(tbl_cnt).powAng = PowPhCr(iN).Ang;
            pTBL(tbl_cnt).powfq_ctrs = PowPhCr(iN).fq_ctrs;
            pTBL(tbl_cnt).powhist_rad = PowPhCr(iN).hist_rad;
            pTBL(tbl_cnt).powAng_p = PowPhCr(iN).Ang_p;
            pTBL(tbl_cnt).powAng_to_shuf_p = PowPhCr(iN).Ang_to_shuf_p;
            pTBL(tbl_cnt).powAng_z = PowPhCr(iN).Ang_z;
            pTBL(tbl_cnt).powsh_hist_rad_mn = PowPhCr(iN).sh_hist_rad_mn;
            pTBL(tbl_cnt).powsh_hist_rad_95ci = PowPhCr(iN).sh_hist_rad_95ci;
            pTBL(tbl_cnt).powsig_ph_locking = PowPhCr(iN).sig_ph_locking;
            pTBL(tbl_cnt).pown_bins_above_shuff = PowPhCr(iN).n_bins_above_shuff;
            
            %                 pTBL(tbl_cnt).Coh_STA_x_sec = PCoh(iN).STA_x_sec  ;
            %                 pTBL(tbl_cnt).Coh_psd_fqs = PCoh(iN).psd_fqs;
            %                 pTBL(tbl_cnt).Coh_circ_rtest_z = PCoh(iN).circ_rtest_z;
            %                 pTBL(tbl_cnt).Coh_circ_rtest_p = PCoh(iN).circ_rtest_p;
            %                 pTBL(tbl_cnt).Coh_circ_rtest_z_center = PCoh(iN).circ_rtest_z_center;
            %                 pTBL(tbl_cnt).Coh_circ_rtest_p_center = PCoh(iN).circ_rtest_p_center;
            %                 pTBL(tbl_cnt).Coh_STA_LFP = PCoh(iN).STA_LFP;
            %                 pTBL(tbl_cnt).Coh_CWT_of_STA = PCoh(iN).CWT_of_STA;
            %                 pTBL(tbl_cnt).Coh_CWT_of_STA_psd = PCoh(iN).CWT_of_STA_psd;
            %                 pTBL(tbl_cnt).Coh_STA_pmtm_psd = PCoh(iN).STA_pmtm_psd;
            %                 pTBL(tbl_cnt).Coh_STA_pmtm_psd_norm = PCoh(iN).STA_pmtm_psd_norm;
            
            %                 pTBL(tbl_cnt).PhPk_STA_x_sec = PhPk(iN).STA_x_sec;
            %                 pTBL(tbl_cnt).PhPk_fq_range = PhPk(iN).fq_range;
            %                 pTBL(tbl_cnt).PhPk_orig_ix = PhPk(iN).orig_ix;
            %                 pTBL(tbl_cnt).PhPk_circ_rtest_z = PhPk(iN).circ_rtest_z;
            %                 pTBL(tbl_cnt).PhPk_circ_rtest_p = PhPk(iN).circ_rtest_p;
            %                 pTBL(tbl_cnt).PhPk_circ_rtest_sh_z = PhPk(iN).circ_rtest_sh_z;
            %                 pTBL(tbl_cnt).PhPk_circ_rtest_sh_p = PhPk(iN).circ_rtest_sh_p;
            %                 pTBL(tbl_cnt).PhPk_circ_z_to_shuff = PhPk(iN).circ_z_to_shuff;
            %                 pTBL(tbl_cnt).PhPk_circ_p_to_shuff = PhPk(iN).circ_p_to_shuff;
            %                 pTBL(tbl_cnt).PhPk_freq_pow_r = PhPk(iN).freq_pow_r;
            %                 pTBL(tbl_cnt).PhPk_freq_pow_p = PhPk(iN).freq_pow_p;
            %                 pTBL(tbl_cnt).PhPk_pk_freq = PhPk(iN).pk_freq;
            %                 pTBL(tbl_cnt).PhPk_pk_pow = PhPk(iN).pk_pow;
            %                 pTBL(tbl_cnt).PhPk_pk_phase = PhPk(iN).pk_phase;
            
            tbl_cnt = tbl_cnt + 1;
        end
        
    end
end

TBL = struct2table(pTBL);

OUT.TBL = TBL;
OUT.pTBL = pTBL;
