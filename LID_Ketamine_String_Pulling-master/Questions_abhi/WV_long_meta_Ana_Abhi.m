%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'WV_long_meta_Abhi';

GP.Analysis_Dir = 'E:\Temp\TempAnaResults\WV_long_meta_Abhi';

PLOT_IT = false;
% adir = fullfile(GP.Analysis_Dir,ses_to_ana);
adir = fullfile(GP.Analysis_Dir);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end
% Initialize
NRN = [];
nrn_cnt = 1;

ALL_WV_long = [];
ALL_WV = [];
ALL_AC_base = [];

for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));

    %get the two SP structs
    SP_WV = Dset.SP;
    SP_WV_long = Dset.SP_WV_long;
  
    % Creating meta table for phase and power modulation information foe
    % each neuron
    for ii = 1:length(SP_WV)
            % general
            NRN(nrn_cnt).NeuronID = SP_WV(ii).fname;
            NRN(nrn_cnt).Rat = Dset.SES.Rat;
            NRN(nrn_cnt).Session = Dset.SES.Session;
            if Dset.SES.RatType == 'Control'
                NRN(nrn_cnt).Group = 'Control';
            elseif Dset.SES.RatType == '6ODHA_LID'
                NRN(nrn_cnt).Group = '6ODHA_LID';
            elseif Dset.SES.RatType == 'SHAM'
                NRN(nrn_cnt).Group = 'Sham';
            end
            NRN(nrn_cnt).Drugs = Dset.SES.Drugs;
            
            NRN(nrn_cnt).Depth_uM = SP_WV(ii).Depth_uM;
            if SP_WV(ii).Depth_uM < 3000
                NRN(nrn_cnt).BrainRegion = 'M1';
            else
                NRN(nrn_cnt).BrainRegion = 'Striatum';
            end
            NRN(nrn_cnt).Hemisphere = SP_WV(ii).Hemisphere;
            NRN(nrn_cnt).Tetrode_Channel = SP_WV(ii).Tetrode;
            
            
            wv = SP_WV(ii).WV.mWV;
            [~,ix] = max(max(wv));
            ALL_WV(nrn_cnt,:) = wv(:,ix)';
            
            wv_long = SP_WV_long(ii).WV_LONG.mn;
            [~,ix_long] = max(max(wv_long'));
            ALL_WV_long(nrn_cnt,:) = wv_long(ix_long,:);
            
            ALL_AC_base(nrn_cnt,:) = Dset.AC_base_post1_post2_post3(ii,:,1);
            
            
            nrn_cnt = nrn_cnt + 1;
    end
end

AC_x_msec = Dset.AC_x_ms;

ALL_WV_long = ALL_WV_long * -1;


% for iu = 1: length(ALL_WV_long)
%     figure(1)
%     clf
%     plot(ALL_WV_long(iu,:))
%     pause
% end

[peak_half_width, peak_to_trough] = Spike_width(ALL_WV_long);
wv_features = [peak_half_width, peak_to_trough];% figure; scatter(peak_half_width, peak_to_trough)
wv_feature_labels = {'half width' 'pk to tr'};

% 32 WV points
WV_point = (1/30000)*1000;
ALL_WV_x_msec = 0:WV_point:WV_point*31;

% Long WV 58 points
ALL_WV_x_msec_long = 0:WV_point:WV_point*57;

plot_it = true;

[N_type,type_lables] = neuron_subtypes(wv_features,wv_feature_labels,ALL_WV_long,ALL_WV_x_msec_long,ALL_AC_base,AC_x_msec,plot_it);


