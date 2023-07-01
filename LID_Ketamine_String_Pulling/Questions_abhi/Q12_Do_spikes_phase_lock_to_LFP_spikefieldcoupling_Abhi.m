function OUT = Q12_Do_spikes_phase_lock_to_LFP_spikefieldcoupling_Abhi()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does LDOPA Create 80-Hz oscillations?
% TODO: Correlate and plot firing rate with phase locking.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DIRS
OUT = [];
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();
%lfp_dir = './LFP';
% lfp_dir = DIRS.LFP_Dir;

DO_TREADMILL = false;
regress_80 = false;
n_ac_bins = 100;
ac_binsize_ms = 4;
ac_x_ms = ((1:n_ac_bins)*ac_binsize_ms) - ac_binsize_ms/2;
PLOT_IT = false;
psd_fqs = 1:.5:250;
fqs_peakcoupling = [35 95];
fq_ctrs = [2 8 22 50 84];
resamp_LFP_freq = 2000;
% plot_type = 'wavelet';
SES = LK_Session_Info();
OUT.SES = SES;
OUT.aborted = false;
OUT.SP = SP;
OUT.n_ac_bins = n_ac_bins;
OUT.ac_binsize_ms = ac_binsize_ms;
tbl_cnt = 1;
OSC = Oscillation_frequencies;
% fq_bands = { 'low_gamma_colgin'  'gamma_80' };

% if exist('LFP_INFO.csv','file')
%     lfp_info = readtable('LFP_INFO.csv');
% else
%     lfp_info = readtable('LFP_INFO_auto_generated.csv');
% end
% lfp_hand_selected = readtable('LFP_information.csv');
% Load the 80 Hz channel and the channel that is good on unlesioned
% hemisphere.
% IX1 = categorical(lfp_hand_selected.signal) == 'HFOReRef1';
% IX2 = categorical(lfp_hand_selected.signal) == 'HFOReRef2';
% fq_bands = {'gamma_50' 'low_gamma_colgin' 'HFO'};

% IX1 = categorical(lfp_hand_selected.signal) == 'Lesioned80HzReRef1';
% IX2 = categorical(lfp_hand_selected.signal) == 'Lesioned80HzReRef2';
% fq_bands = {'gamma_50' 'low_gamma_colgin' 'gamma_80'};
% fq_bands = {'beta' 'low_gamma_colgin' 'gamma_80'};
% fq_bands = {'gamma_80'};
fq_bands = {'delta' 'theta' 'beta' 'gamma_50' 'gamma_80'};

%   IX1 = categorical(lfp_hand_selected.signal) == 'LesionedBetaReRef1';
%   IX2 = categorical(lfp_hand_selected.signal) == 'LesionedBetaReRef2';
%   fq_bands = {'beta3' 'beta' 'theta' };

%  IX1 = categorical(lfp_hand_selected.signal) == 'UnlesionedKetamine50HzReRef1';
%  IX2 = categorical(lfp_hand_selected.signal) == 'UnlesionedKetamine50HzReRef2';
% sigfile = lfp_hand_selected.filename(IX1);
% reffile = lfp_hand_selected.filename(IX2);
% sigfile = lfp_info.fname(lfp_info.good_sig==1);
% [LFP] = LK_Load_and_Clean_LFP(sigfile{1}, lfp_info.fname(lfp_info.good_ref==1));
% [LFP] = LK_Load_and_Clean_LFP(lfp_dir, sigfile{1}, reffile);
% load('./Processed_Data/best_channels_gamma_80.mat','OUT') %Stephen's code
% Load LFP file 
LFPfiles = find_files(fullfile('M1*.mat'));
LFP = LK_Load_and_Clean_LFP_Abhi('',LFPfiles{1});
%% Load the best non-reref.
% if exist(fullfile('LFP',OUT.best_non_reref),'file')
%     % load if stored locally.
%     LFP = LK_Load_and_Clean_LFP_Abhi('./LFP',OUT.best_non_reref);
% else
%     global DIRS
%     lfp_dir = fullfile(DIRS.LFP_Dir,SES.rat_str,SES.session_str,'LFP');
%     LFP = LK_Load_and_Clean_LFP(lfp_dir, OUT.best_non_reref);
% end

% IMU = LK_Load_and_Process_IMU;
% POS = LK_Load_and_Clean_POS;
% determine what events are available for alignment.
te = []; cnt = 1;
intervals_around_drug_min = [];
if any(E.EventID == 'KetInjectionStart')
    te{cnt} = 'KET';
    e_start_end_uS(cnt,1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
    e_start_end_uS(cnt,2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
    intervals_around_drug_min{cnt} = [-25 -5; 2 30; 60 80];
    
    cnt = cnt + 1;
end
if any(E.EventID == 'LDOPAInjectionStart')
    te{cnt} = 'LDO';
    e_start_end_uS(cnt,1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
    e_start_end_uS(cnt,2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
%     intervals_around_drug_min{cnt} = [-25 -5; 25 55; 65 90];
    intervals_around_drug_min{cnt} = [-25 -5; 5 15; 25 55; 65 90]; % added 9/26/22
    
    cnt = cnt + 1;
end
if DO_TREADMILL
    if any(E.EventID == 'TreadmillStart')
        te{cnt} = 'Tread';
        tmp1 = E.MinFromStart(E.EventID == 'TreadmillStart');
        %     tmp2 = E.MinFromStart(E.EventID == 'TreadmillStop');
        %
        e_start_end_uS(cnt,1) = tmp1(1)*60*1e6;
        e_start_end_uS(cnt,2) = tmp1(1)*60*1e6 + 10e6;
        
        %     e_start_end_uS(cnt,1) = E.MinFromStart(E.EventID == 'TreadmillStart')*60*1e6 +1;
        %     e_start_end_uS(cnt,2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        intervals_around_drug_min{cnt} = [-20 -5; 2 22; 25 40];
        
        cnt = cnt + 1;
    end
end
% Create some filters.
F = SPEC_create_filters(fq_bands,LFP.sFreq);
% For testing, make an artificial neuron phase locked to a frequency.
if 0
    SP(end+1) = SP(end);
    TS{end+1} = [];
    ix1 = binsearch(LFP.t_uS , e_start_end_uS(end,2));
    ix2 = binsearch(LFP.t_uS , e_start_end_uS(end,2)+60e6);
    
    L = filtfilt(F_adj{1},LFP.LFP(ix1:ix2))';
    t = LFP.t_uS(ix1:ix2);
    [pk,ix] = findpeaks(L);
    
    SP(end).t_uS = t(ix(pk>prctile(pk,60)));
    TS{end} = SP(end).t_uS;
end

pTBL = struct();

% Filter and upsample LFP
LFP_upsamp_uS = LFP.t_uS(1):1e6/resamp_LFP_freq:LFP.t_uS(end);
LFP_upsamp = [];
for iFreq = 1:length(fq_bands)
    % Filter data bandpass for each freq and see if spikes phase
    % lock.
    
    LFP_filt = filtfilt(F{iFreq},LFP.LFP);
    LFP_upsamp(:,iFreq) = interp1(LFP.t_uS,LFP_filt,LFP_upsamp_uS,'spline');
    
end

%% Creating envelopes of the adjacent bands to gamma 80 to get a measure 
% of focal 80 but subtracting the adjacent bands from 80 later in the code

% adj_bands = {[68 72] [98 102]};

% F_adj = SPEC_create_filters(adj_bands,LFP.sFreq);
% for b = 1:length(adj_bands)
% 
%     L_adj(:,b) = filtfilt(F_adj{b},LFP.LFP);
% %     % power envelope
% %     ss = abs(hilbert(L_adj(:,b)));
% %     ss = ss - mean(ss);
% %     Lp_adj(:,b) = ss;
%     
% end
for iTE = 1:length(te)
    
    TIMES.(te{iTE}).EventStartEndUsec(1) = e_start_end_uS(iTE,1);
    TIMES.(te{iTE}).EventStartEndUsec(2) = e_start_end_uS(iTE,2);
    TIMES.(te{iTE}).IntervalUsec(1,1) = TIMES.(te{iTE}).EventStartEndUsec(1) + intervals_around_drug_min{iTE}(1,1)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(1,2) = TIMES.(te{iTE}).EventStartEndUsec(1) + intervals_around_drug_min{iTE}(1,2)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(2,1) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(2,1)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(2,2) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(2,2)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(3,1) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(3,1)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(3,2) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(3,2)*60*1e6;
    
    if length(intervals_around_drug_min{iTE}) > 3
        TIMES.(te{iTE}).IntervalUsec(4,1) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(4,1)*60*1e6;
        TIMES.(te{iTE}).IntervalUsec(4,2) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(4,2)*60*1e6;
    else
    end
    
    TIMES.(te{iTE}).PeriEventUsec = [intervals_around_drug_min{iTE}(1,1) intervals_around_drug_min{iTE}(end,end) ];
    
    OUT.TM.(te{iTE}) = TIMES;
    
    
        for ii = 1:Rows(TIMES.(te{iTE}).IntervalUsec)
        ix1 = binsearch(LFP_upsamp_uS , TIMES.(te{iTE}).IntervalUsec(ii,1));
        ix2 = binsearch(LFP_upsamp_uS , TIMES.(te{iTE}).IntervalUsec(ii,2));
        L = [];   
        L = [LFP_upsamp_uS(ix1:ix2)' LFP_upsamp(ix1:ix2,:)];
            
%             L(:,1) = LFP_upsamp_uS(ix1:ix2);
% 
%             for iFreq = 1:length(fq_bands)
%                 % Filter data bandpass for each freq and see if spikes phase
%                 % lock.
%               
%                 L(:,iFreq+1) = LFP_upsamp(ix1:ix2,iFreq);
% %                 
%             end
%                 
        TSr = Restrict(TS,LFP_upsamp_uS(ix1),LFP_upsamp_uS(ix2));
        
%         Ang =  angle(hilbert(L(:,4)));

% % 
%         [PowPhCr]  = SPEC_spike_field_coupling(TS,[LFP.t_uS(ix1:ix2), Lp],[8 22 50 84 150],'thresh_prctile_of_power',80);
% 
        [PhC,phinfo] = SPEC_spike_field_coupling(TSr,L,fq_ctrs,'thresh_prctile_of_power',[80 10]);
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

            for iN = 1:length(SP)
                % meta stuff
                pTBL(tbl_cnt).NeuronID = SP(iN).fname;
%                 pTBL(tbl_cnt).UNeuronID = SES.rat + 1000*SES.session + 100000*iN;
                pTBL(tbl_cnt).Session = SES.session;
                pTBL(tbl_cnt).Rat = SES.rat;
                pTBL(tbl_cnt).Condition = te{iTE};
%                 pTBL(tbl_cnt).FreqBand = fq_bands{iFreq};
                pTBL(tbl_cnt).Interval = ii;
                pTBL(tbl_cnt).Tetrode = SP(iN).Tetrode;
                pTBL(tbl_cnt).Hemisphere = SP(iN).Hemisphere;
                pTBL(tbl_cnt).Depth_uM = SP(iN).Depth_uM;
                pTBL(tbl_cnt).WV = single(SP(iN).WV.mWV);

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
                pTBL(tbl_cnt).TS = TS(iN,1);
                
                % Calculated stuff
% %                 pTBL(tbl_cnt).Ang = PhC(iN).Ang;
                pTBL(tbl_cnt).fq_ctrs = PhC(iN).fq_ctrs;
                pTBL(tbl_cnt).hist_rad = PhC(iN).hist_rad;
                % pTBL(tbl_cnt).Ang_ts = PhC(iN).Ang_ts;
                pTBL(tbl_cnt).Ang_p = PhC(iN).Ang_p;
                pTBL(tbl_cnt).Ang_to_shuf_p = PhC(iN).Ang_to_shuf_p;
                pTBL(tbl_cnt).Ang_to_shuf_F = PhC(iN).Ang_to_shuf_F;
                pTBL(tbl_cnt).Ang_z = PhC(iN).Ang_z;
                pTBL(tbl_cnt).sh_hist_rad_mn = PhC(iN).sh_hist_rad_mn;
                pTBL(tbl_cnt).sh_hist_rad_95ci = PhC(iN).sh_hist_rad_95ci;
                pTBL(tbl_cnt).sig_ph_locking = PhC(iN).sig_ph_locking;
                pTBL(tbl_cnt).n_bins_above_shuff = PhC(iN).n_bins_above_shuff;
                pTBL(tbl_cnt).good_pow_intervals = PhC(iN).good_pow_intervals;
                
                
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
% %                 pTBL(tbl_cnt).powAng = PowPhCr(iN).Ang;
%                 pTBL(tbl_cnt).powfq_ctrs = PowPhCr(iN).fq_ctrs;
%                 pTBL(tbl_cnt).powhist_rad = PowPhCr(iN).hist_rad;
%                 pTBL(tbl_cnt).powAng_p = PowPhCr(iN).Ang_p;
%                 pTBL(tbl_cnt).powAng_to_shuf_p = PowPhCr(iN).Ang_to_shuf_p;
%                 pTBL(tbl_cnt).powAng_z = PowPhCr(iN).Ang_z;
%                 pTBL(tbl_cnt).powsh_hist_rad_mn = PowPhCr(iN).sh_hist_rad_mn;
%                 pTBL(tbl_cnt).powsh_hist_rad_95ci = PowPhCr(iN).sh_hist_rad_95ci;
%                 pTBL(tbl_cnt).powsig_ph_locking = PowPhCr(iN).sig_ph_locking;
%                 pTBL(tbl_cnt).pown_bins_above_shuff = PowPhCr(iN).n_bins_above_shuff;
%                 
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
% TBL.Condition = categorical(TBL.Condition);
% TBL.FreqBand = categorical(TBL.FreqBand);
% TBL.Hemisphere = categorical(TBL.Hemisphere);
% Plot
% Spec
if 0
    for ii = 1:3
        figure(1010)
        subplot(3,1,ii)
        ix1 = binsearch(LFP.t_uS , TIMES.(te{iTE}).IntervalUsec(ii,1));
        ix2 = binsearch(LFP.t_uS , TIMES.(te{iTE}).IntervalUsec(ii,2));
        spectrogram(LFP.LFP(ix1:ix2),LFP.sFreq*4,LFP.sFreq*2,256,LFP.sFreq,'yaxis')
        figure(1011)
        subplot(3,1,ii)
        e = envelope_cowen(filtfilt(F_adj{1},LFP.LFP(ix1:ix2)));
        plot(e);
        m(ii) = nanmean(e);
        
        figure(1013)
        subplot(3,1,ii)
        pwelch(LFP.LFP(ix1:ix2),LFP.sFreq*4,LFP.sFreq*2,256,LFP.sFreq)
        
    end
    
    figure(1011)
    equalize_y_axes
    figure
    bar(m)
    % equalize_color_axes
end
%%


if PLOT_IT
    clrs = lines(3);
    % find neurons with a selective response
    for iTE = 1:length(te)
        for iFq = 1:length(fq_bands)
            CIX = TBL.Condition == te{iTE} & TBL.FreqBand == fq_bands{iFq};
            GOOD_NEURONS  = unique(TBL.NeuronID(TBL.Ang_p < 0.05) );
            ix  = find(ismember(TBL.NeuronID,GOOD_NEURONS) & CIX);
            for jj = 1:length(ix)
                row_ix = ix(jj);
                iN = TBL.NeuronID(ix(jj));
                ii = TBL.Interval(row_ix);
                fq = char(TBL.FreqBand(row_ix));
                
                str = sprintf('%d %s %s',iN,fq,TBL.Condition(row_ix) );
                
                figure(iN + 100*iTE + 1000*iFq)
                if ~isnan(TBL.AC{row_ix})
                    v1 = TBL.sh_hist_rad_mn{row_ix};
                    v2 = TBL.powsh_hist_rad_mn{row_ix};
                    if ~isempty(v1)
                        
                        subplot(3,3,ii)
                        polarhistogram('BinEdges',phinfo.rad_edges,'BinCounts',v1 ,'DisplayStyle','stairs','EdgeColor','k')
                        hold on
                        polarhistogram('BinEdges',phinfo.rad_edges,'BinCounts',TBL.hist_rad{row_ix} ,'FaceAlpha',.2,'EdgeColor',clrs(ii,:),'FaceColor',clrs(ii,:))
                        
                        title(sprintf ('ph %0.3f', TBL.Ang_p(row_ix)))                        
                    end
                    if ~isempty(v2)
                        subplot(3,3,ii+3)
                        polarhistogram('BinEdges',phinfo.rad_edges,'BinCounts',v2 ,'DisplayStyle','stairs','EdgeColor','k')
                        hold on
                        polarhistogram('BinEdges',phinfo.rad_edges,'BinCounts',TBL.powhist_rad{row_ix},'FaceAlpha',.2,'EdgeColor',clrs(ii,:),'FaceColor',clrs(ii,:))
                        
                        title(sprintf ('pw %0.3f', TBL.powAng_p(row_ix)))
                    end
                    subplot(3,3,ii+6)
                    plot(ac_x_ms,TBL.AC{row_ix},'Color',clrs(ii,:))
                    axis tight
                    pubify_figure_axis
                    xlabel('ms')
                end
                %         TBL.NeuronID
                sgtitle(str)
            end
        end
    end
end

%% PLot autocorrs.
%
%         for iTE = 1:length(te)
%             figure
%             for iInt = 1:3
%                 GIX = ~isnan(squeeze(AC.(te{iTE})(1,iInt,:)));
%                 subplot(2,3,iInt)
%                 imagesc(ac_x_ms,[],squeeze(AC.(te{iTE})(:,iInt,GIX))')
%                 title(te{iTE})
%                 subplot(2,3,iInt+3)
%                 plot_confidence_intervals(ac_x_ms,squeeze(AC.(te{iTE})(:,iInt,GIX))')
%             end
%         end
%     end
% OUT.AC = AC;
OUT.SPmod = SP;
OUT.TBL = TBL;
OUT.pTBL = pTBL;


     
            
            
