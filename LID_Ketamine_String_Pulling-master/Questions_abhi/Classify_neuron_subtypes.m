%% Analyze data across sessions.
% For 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
ses_to_ana = 'Get_AutoCorr';
Analysis_Dir = 'E:\Temp\TempAnaResults\Q1_Sal_Ket';

adir = fullfile(Analysis_Dir,ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end

ALL_AC= [];
ALL_WV= [];
AC_x_msec = [];
NRN = [];
nrn_cnt = 1;

for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));
    
     for ii = 1:Rows(Dset.AC)
        ALL_AC(nrn_cnt,:) = Dset.AC(ii,:);
        
        wv = Dset.SP(ii).WV.mWV;
        [~,ix] = max(max(wv));
        ALL_WV(nrn_cnt,:) = wv(:,ix)';
        
        NRN(nrn_cnt).NeuronID = Dset.SP(ii).fname;
        NRN(nrn_cnt).Rat = Dset.SES.Rat;
        NRN(nrn_cnt).Session = Dset.SES.Session;
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
        if ~isempty(Dset.SP(ii).Tetrode_Channels)
            NRN(nrn_cnt).Tetrode_Channel_1 = Dset.SP(ii).Tetrode_Channels(1);
        else
            NRN(nrn_cnt).Tetrode_Channel_1 = nan;
        end
        
        nrn_cnt = nrn_cnt + 1;
    end
    
end

TBL = struct2table(NRN);

IX_M1 = categorical(TBL.BrainRegion) == 'M1';

AC_x_msec = Dset.AC_x_ms; 
WV_point = (1/30000)*1000;
ALL_WV_x_msec = 0:WV_point:WV_point*31;

% Get the WV features
[peak_half_width, peak_to_trough] = Spike_width(ALL_WV);

wv_features = [peak_half_width, peak_to_trough];%
wv_feature_labels = {'half width' 'pk to tr'};

plot_it = true;

[IX_M1_type,type_lables] = neuron_subtypes(wv_features(IX_M1,:),wv_feature_labels,ALL_WV(IX_M1,:),ALL_WV_x_msec,ALL_AC(IX_M1,:),AC_x_msec,plot_it);
 
g = kmeans(nan_to_val(wv_features),3); % c
figure;gscatter(wv_features(:,1),wv_features(:,2),g)
% tsne plot with wv_features and AC

