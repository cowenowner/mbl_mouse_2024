%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q12_Do_spikes_phase_lock_to_LFP_spikefieldcoupling_Ana';
ses_to_ana = 'Q12_Do_spikes_phase_lock_to_LFP_spikefieldcoupling_Abhi';

GP.Analysis_Dir = 'E:\Temp_Abhi';

PLOT_IT = true;
adir = fullfile(GP.Analysis_Dir,ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end
cm = lines(5);

% Initialize
NRN = [];
nrn_cnt = 1;
wv_cnt = 1;
ALL_WV = [];
ALL_AC = [];
SP_WV = [];
for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));

    % Convert the TBL in Dset file to a structure
    sTBL = table2struct(Dset.TBL);
    
    for ip = 1:length(Dset.SP)
        wv = Dset.SP(ip).WV_LONG.mn';
        [~,ix] = max(max(wv));
        ALL_WV(wv_cnt,:) = wv(:,ix)';
        
        SP_WV(wv_cnt).NeuronID = Dset.SP(ip).fname;
        SP_WV(wv_cnt).Rat = Dset.SES.Rat;
        SP_WV(wv_cnt).Session = Dset.SES.Session;
        SP_WV(wv_cnt).Drugs = Dset.SES.Drugs;
        if Dset.SES.RatType == '6ODHA_LID'
            SP_WV(wv_cnt).Group = '6OHDA_LID';
        elseif Dset.SES.RatType == 'SHAM'
            SP_WV(wv_cnt).Group = 'SHAM';
        end
        if Dset.SP(ip).Depth_uM < 3000
            SP_WV(wv_cnt).BrainRegion = 'M1';
        else
            SP_WV(wv_cnt).BrainRegion = 'Striatum';
        end
        wv_cnt = wv_cnt + 1;
    end
    
    % Creating meta table for phase and power modulation information foe
    % each neuron
    for ii = 1:length(sTBL)
            % general
            NRN(nrn_cnt).NeuronID = sTBL(ii).NeuronID;
            NRN(nrn_cnt).Rat = Dset.SES.Rat;
            NRN(nrn_cnt).Session = Dset.SES.Session;
            if Dset.SES.RatType == '6ODHA_LID'
                NRN(nrn_cnt).Group = '6OHDA_LID';
            elseif Dset.SES.RatType == 'SHAM'
                NRN(nrn_cnt).Group = 'SHAM';
            end
            NRN(nrn_cnt).Drugs = Dset.SES.Drugs;
            
            NRN(nrn_cnt).Depth_uM = sTBL(ii).Depth_uM;
            if sTBL(ii).Depth_uM < 3000
                NRN(nrn_cnt).BrainRegion = 'M1';
            else
                NRN(nrn_cnt).BrainRegion = 'Striatum';
            end
            NRN(nrn_cnt).Hemisphere = sTBL(ii).Hemisphere;
            NRN(nrn_cnt).Tetrode_Channel = sTBL(ii).Tetrode;
%             NRN(nrn_cnt).WV = sTBL(ii).WV;
            NRN(nrn_cnt).AC = sTBL(ii).AC;
            NRN(nrn_cnt).Condition = sTBL(ii).Condition;
            NRN(nrn_cnt).Interval = sTBL(ii).Interval;
            
            NRN(nrn_cnt).fq_ctrs = sTBL(ii).fq_ctrs;
            NRN(nrn_cnt).hist_rad = sTBL(ii).hist_rad;
            NRN(nrn_cnt).Ang_p = sTBL(ii).Ang_p;
            NRN(nrn_cnt).Ang_to_shuf_p = sTBL(ii).Ang_to_shuf_p;
            NRN(nrn_cnt).Ang_z = sTBL(ii).Ang_z;
            NRN(nrn_cnt).sh_hist_rad_mn = sTBL(ii).sh_hist_rad_mn;
            NRN(nrn_cnt).sh_hist_rad_95ci = sTBL(ii).sh_hist_rad_95ci;
            NRN(nrn_cnt).sig_ph_locking = sTBL(ii).sig_ph_locking;
            NRN(nrn_cnt).n_bins_above_shuff = sTBL(ii).n_bins_above_shuff;
%             NRN(nrn_cnt).LFP_filt = sTBL(ii).LFP_uS_filt;
            NRN(nrn_cnt).good_pow_intervals = sTBL(ii).good_pow_intervals;
            NRN(nrn_cnt).TS = sTBL(ii).TS;
            %         NRN(nrn_cnt).powhist_rad = sTBL(ii).powhist_rad;
            %         NRN(nrn_cnt).powAng_p = sTBL(ii).powAng_p;
            %         NRN(nrn_cnt).powAng_to_shuf_p = sTBL(ii).powAng_to_shuf_p;
            %         NRN(nrn_cnt).powAng_z = sTBL(ii).powAng_z;
            %         NRN(nrn_cnt).powsh_hist_rad_mn = sTBL(ii).powsh_hist_rad_mn;
            %         NRN(nrn_cnt).powsh_hist_rad_95ci = sTBL(ii).powsh_hist_rad_95ci;
            %         NRN(nrn_cnt).powsig_ph_locking = sTBL(ii).powsig_ph_locking;
            %         NRN(nrn_cnt).pown_bins_above_shuff = sTBL(ii).pown_bins_above_shuff;
            %
            %         NRN(nrn_cnt).Coh_STA_x_sec = sTBL(ii).Coh_STA_x_sec  ;
            %         NRN(nrn_cnt).Coh_psd_fqs = sTBL(ii).Coh_psd_fqs;
            %         NRN(nrn_cnt).Coh_circ_rtest_z = sTBL(ii).Coh_circ_rtest_z;
            %         NRN(nrn_cnt).Coh_circ_rtest_p = sTBL(ii).Coh_circ_rtest_p;
            %         NRN(nrn_cnt).Coh_circ_rtest_z_center = sTBL(ii).Coh_circ_rtest_z_center;
            %         NRN(nrn_cnt).Coh_circ_rtest_p_center = sTBL(ii).Coh_circ_rtest_p_center;
            %         NRN(nrn_cnt).Coh_STA_LFP = sTBL(ii).Coh_STA_LFP;
            %         NRN(nrn_cnt).Coh_CWT_of_STA = sTBL(ii).Coh_CWT_of_STA;
            %         NRN(nrn_cnt).Coh_CWT_of_STA_psd = sTBL(ii).Coh_CWT_of_STA_psd;
            %         NRN(nrn_cnt).Coh_STA_pmtm_psd = sTBL(ii).Coh_STA_pmtm_psd;
            %         NRN(nrn_cnt).Coh_STA_pmtm_psd_norm = sTBL(ii).Coh_STA_pmtm_psd_norm;
            %
            
            %         NRN(nrn_cnt).PhPk_STA_x_sec = sTBL(ii).PhPk_STA_x_sec;
            %         NRN(nrn_cnt).PhPk_fq_range = sTBL(ii).PhPk_fq_range;
            %         NRN(nrn_cnt).PhPk_circ_rtest_z = sTBL(ii).PhPk_circ_rtest_z;
            %         NRN(nrn_cnt).PhPk_circ_rtest_p = sTBL(ii).PhPk_circ_rtest_p;
            %         NRN(nrn_cnt).PhPk_circ_rtest_sh_z = sTBL(ii).PhPk_circ_rtest_sh_z;
            %         NRN(nrn_cnt).PhPk_circ_rtest_sh_p = sTBL(ii).PhPk_circ_rtest_sh_p;
            %         NRN(nrn_cnt).PhPk_circ_z_to_shuff = sTBL(ii).PhPk_circ_z_to_shuff;
            %         NRN(nrn_cnt).PhPk_circ_p_to_shuff = sTBL(ii).PhPk_circ_p_to_shuff;
            %         NRN(nrn_cnt).PhPk_freq_pow_r = sTBL(ii).PhPk_freq_pow_r;
            %         NRN(nrn_cnt).PhPk_freq_pow_p = sTBL(ii).PhPk_freq_pow_p;
            %         NRN(nrn_cnt).PhPk_pk_freq = sTBL(ii).PhPk_pk_freq;
            %         NRN(nrn_cnt).PhPk_pk_pow = sTBL(ii).PhPk_pk_pow;
            %         NRN(nrn_cnt).PhPk_pk_phase = sTBL(ii).PhPk_pk_phase;
            
            ALL_AC(nrn_cnt,:) = sTBL(ii).AC;
            
            nrn_cnt = nrn_cnt + 1;

    end
end
% save to table
TBL = struct2table(NRN);
for ih = 1:length(ALL_WV)
    SP_WV(ih).WV = ALL_WV(ih,:);
    
end
SP_TBL = struct2table(SP_WV);

% did a re classification of good spikes based on waveforms and got rid of
% weird looking waveforms
SP_TBL = SP_TBL(~IX_bad_spikes_all,:);
ALL_WV = SP_TBL.WV;


IX_M1 = categorical(TBL.BrainRegion) == 'M1';
TBL_M1 = TBL(IX_M1,:);
TBL_M1 = TBL_M1(~IX_bad_spikes_all,:);

% Neuron type analysis
%
IX_SP_M1 = categorical(SP_TBL.BrainRegion) == 'M1';
WV_base_M1 = ALL_WV(IX_SP_M1,:);
WV_base_M1 = WV_base_M1 * -1;

IX_Ket1_M1 = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 1 & categorical(TBL_M1.BrainRegion) == 'M1'; 
IX_LK_LS = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' | categorical(TBL_M1.Drugs) == 'LDOPA&Saline';
IX_LDO1_M1 = IX_LK_LS & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 1 & categorical(TBL_M1.BrainRegion) == 'M1';
IX_WV_analysis = IX_Ket1_M1 | IX_LDO1_M1;
TBL_WV = TBL_M1(IX_WV_analysis,:);

ALL_AC_M1 = ALL_AC(IX_M1,:);
AC_base_M1 = ALL_AC_M1(IX_WV_analysis,:);
% 
% for iw = 1:size(TBL_WV)
%     
%     wv = TBL_WV.WV{iw};
%     [~,ix] = max(max(wv));
%     ALL_WV(iw,:) = wv(:,ix)';
% %     AC =  TBL_WV.AC{iw};
%     ALL_AC_base(iw,:) = TBL_WV.AC{iw};
% end

n_ac_bins = 100;
ac_binsize_ms = 4;
AC_x_msec = ((1:n_ac_bins)*ac_binsize_ms) - ac_binsize_ms/2;

WV_point = (1/30000)*1000;
ALL_WV_x_msec = 0:WV_point:WV_point*57;
% creating waveform features
[peak_half_width, peak_to_trough] = Spike_width(WV_base_M1);
wv_features = [peak_half_width, peak_to_trough];% figure; scatter(peak_half_width, peak_to_trough)
wv_feature_labels = {'half width' 'pk to tr'};
plot_it = true;

figure; scatter(peak_half_width, peak_to_trough)

[N_type,type_lables] = neuron_subtypes(wv_features,wv_feature_labels,WV_base_M1,ALL_WV_x_msec,AC_base_M1,AC_x_msec,plot_it);
TBL_WV.Neuron_type = N_type;

for iy = 1:size(TBL_M1)
       ID = TBL_M1.NeuronID{iy};
       d = find(categorical(TBL_WV.NeuronID) == ID);
       nt = TBL_WV.Neuron_type(d);
       TBL_M1.Neuron_type(iy) = nt;
end
for ii = 1:length(ALL_WV)
figure
plot(ALL_WV(ii,:))
end