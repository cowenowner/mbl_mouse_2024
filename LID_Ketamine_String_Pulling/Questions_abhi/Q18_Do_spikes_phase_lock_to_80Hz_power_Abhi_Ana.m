%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q18_Do_spikes_phase_lock_to_80Hz_power_Abhi_Ana';
ses_to_ana = 'Q18_Do_spikes_phase_lock_to_80Hz_power_Abhi';

GP.Analysis_Dir = 'H:\Temp';

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

    % Get the structure
    pTBL = Dset.pTBL;
    
    % Creating meta table for phase and power modulation information foe
    % each neuron
    if Dset.aborted == false
        for ii = 1:length(pTBL)
            % general
            NRN(nrn_cnt).NeuronID = pTBL(ii).NeuronID;
            NRN(nrn_cnt).Rat = Dset.SES.Rat;
            NRN(nrn_cnt).Session = Dset.SES.Session;
            if Dset.SES.RatType == '6ODHA_LID'
                NRN(nrn_cnt).Group = '6OHDA_LID';
            elseif Dset.SES.RatType == 'SHAM'
                NRN(nrn_cnt).Group = 'SHAM';
            end
            NRN(nrn_cnt).Drugs = Dset.SES.Drugs;
            
            NRN(nrn_cnt).Depth_uM = pTBL(ii).Depth_uM;
            if pTBL(ii).Depth_uM < 3000
                NRN(nrn_cnt).BrainRegion = 'M1';
            else
                NRN(nrn_cnt).BrainRegion = 'Striatum';
            end
            NRN(nrn_cnt).Hemisphere = pTBL(ii).Hemisphere;
            NRN(nrn_cnt).Tetrode_Channel = pTBL(ii).Tetrode;
            %             NRN(nrn_cnt).WV = sTBL(ii).WV;            
            NRN(nrn_cnt).Interval = pTBL(ii).Interval;
            NRN(nrn_cnt).LFP_filename = pTBL(ii).LFP_filename;
            NRN(nrn_cnt).AC = pTBL(ii).AC;
            
            NRN(nrn_cnt).fq_ctrs = pTBL(ii).powfq_ctrs;
            NRN(nrn_cnt).hist_rad = pTBL(ii).powhist_rad;
            NRN(nrn_cnt).Ang_p = pTBL(ii).powAng_p;
            NRN(nrn_cnt).Ang_to_shuf_p = pTBL(ii).powAng_to_shuf_p;
            NRN(nrn_cnt).Ang_z = pTBL(ii).powAng_z;
            NRN(nrn_cnt).sh_hist_rad_mn = pTBL(ii).powsh_hist_rad_mn;
            NRN(nrn_cnt).sh_hist_rad_95ci = pTBL(ii).powsh_hist_rad_95ci;
            NRN(nrn_cnt).sig_ph_locking = pTBL(ii).powsig_ph_locking;
            NRN(nrn_cnt).n_bins_above_shuff = pTBL(ii).pown_bins_above_shuff;
            %             NRN(nrn_cnt).LFP_filt = sTBL(ii).LFP_uS_filt;
            NRN(nrn_cnt).TS = pTBL(ii).TSr;
            
            
            ALL_AC(nrn_cnt,:) = pTBL(ii).AC;
            wv = Dset.SP(ip).WV_LONG.mn';
            [~,ix] = max(max(wv));
            ALL_WV(nrn_cnt,:) = wv(:,ix)';
            
            nrn_cnt = nrn_cnt + 1;
            
        end
    else
    end
end
% save to table
TBL = struct2table(NRN);