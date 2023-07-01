function OUT = Q14_Do_neurons_increase_firing_at_onset_of_80Hzburst_Abhi()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abhi 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OUT = [];

[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;

binsize_TS_ms = 20; % 20 ms bins for FR 
n_ac_bins = 200;
ac_binsize_ms = 2;
ac_x_ms = ((1:n_ac_bins)*ac_binsize_ms) - ac_binsize_ms/2;
PLOT_IT = false;
psd_fqs = 1:.5:250;
fqs_peakcoupling = [35 95];
fq_ctrs = 84;
% plot_type = 'wavelet';
SES = LK_Session_Info();
OUT.SES = SES;
OUT.aborted = false;
OUT.n_ac_bins = n_ac_bins;
OUT.ac_binsize_ms = ac_binsize_ms;
OUT.binsize_TS_ms = binsize_TS_ms;
tbl_cnt = 1;

min_spikes = 20; 
thresh_prctile_of_power = [80 10];

minimum_dur_sec = [];
if ~isempty(fq_ctrs)
    for iC = 1:length(fq_ctrs)
       minimum_dur_sec(iC) = (1/fq_ctrs(iC))*thresh_prctile_of_power(2);
    end
end

minimum_dur_in_samples = 200;

fq_bands = {'gamma_80'};
LFPfiles = find_files(fullfile('M1*.mat'));
LFP = LK_Load_and_Clean_LFP_Abhi('',LFPfiles{1});
LFP_sfreq = LFP.sFreq;

% determine what events are available for alignment.
% te = []; cnt = 1;
intervals_around_drug_min = [];

if any(E.EventID == 'LDOPAInjectionStart')
%     te{cnt} = 'LDO';
    e_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
    e_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
%     intervals_around_drug_min{cnt} = [-25 -5; 25 55; 65 90];
    intervals_around_drug_min = [-25 -5; 5 15; 25 55]; % added 9/26/22
    
%     cnt = cnt + 1;
end

% Create 80Hz filter.
F = SPEC_create_filters(fq_bands,LFP.sFreq);

bTBL = struct();

%%
TIMES.EventStartEndUsec(1) = e_start_end_uS(1);
TIMES.EventStartEndUsec(2) = e_start_end_uS(2);
TIMES.IntervalUsec(1,1) = TIMES.EventStartEndUsec(1) + intervals_around_drug_min(1,1)*60*1e6;
TIMES.IntervalUsec(1,2) = TIMES.EventStartEndUsec(1) + intervals_around_drug_min(1,2)*60*1e6;
TIMES.IntervalUsec(2,1) = TIMES.EventStartEndUsec(2) + intervals_around_drug_min(2,1)*60*1e6;
TIMES.IntervalUsec(2,2) = TIMES.EventStartEndUsec(2) + intervals_around_drug_min(2,2)*60*1e6;
TIMES.IntervalUsec(3,1) = TIMES.EventStartEndUsec(2) + intervals_around_drug_min(3,1)*60*1e6;
TIMES.IntervalUsec(3,2) = TIMES.EventStartEndUsec(2) + intervals_around_drug_min(3,2)*60*1e6;

TIMES.PeriEventUsec = [intervals_around_drug_min(1,1) intervals_around_drug_min(end,end) ];

OUT.TM = TIMES;

for ii = 1:Rows(TIMES.IntervalUsec)
    ix1 = binsearch(LFP.t_uS , TIMES.IntervalUsec(ii,1));
    ix2 = binsearch(LFP.t_uS , TIMES.IntervalUsec(ii,2));
    L = [];
    %             Lp = [];
    %             Lu = [];
    for iFreq = 1:length(fq_bands)
        % Filter data bandpass for each freq and see if spikes phase
        % lock.
        L(:,1) = LFP.t_uS(ix1:ix2);
        L(:,iFreq+1) = filtfilt(F{iFreq},LFP.LFP(ix1:ix2));
        %

    end
%     edges = [];
%     TS_bin = [];
%     bins_uS = [];
%     bin_ctrs_TS = [];
    
    TSr = Restrict(TS,LFP.t_uS(ix1),LFP.t_uS(ix2));
    
    edges = TIMES.IntervalUsec(ii,1):binsize_TS_ms*1000:TIMES.IntervalUsec(ii,2);
    [TS_bin,bins_uS] = Bin_ts_array(TS, edges);
    
    bin_ctrs_TS = (bins_uS(:,1) + bins_uS(:,2))/2;
    
    interval_s = abs((TIMES.IntervalUsec(ii,1)-TIMES.IntervalUsec(ii,2))/1e6);
    TS_sfreq = length(bin_ctrs_TS)/interval_s;
    
    % get power intervals greater than a certain percentile in the oscillations 
    good_pow_intervals = cell(Cols(L)-1,1);
    for iC = 2:Cols(L)
        if ~isempty(thresh_prctile_of_power)
            PW = envelope_cowen(abs(L(:,iC)));
            th = prctile(PW(1:4:end),thresh_prctile_of_power(1));
            I = convn(PW>th,ones(1,minimum_dur_in_samples),'same')>0;
            %         good_pow_intervals{iC-1} = find_intervals([DATA(:,1) I],0);
            pow_intervals = find_intervals([L(:,1) I],0);
            GIX = pow_intervals(:,2)-pow_intervals(:,1) >= minimum_dur_sec(iC-1)*1e6;
            good_pow_intervals{iC-1} = pow_intervals(GIX,:);
            % make sure each LFP interval is longer than some minimum duration
            % so that there is enough data to do a meaninful analysis.
            %         GIX = I(:,2)-I(:,1) >= minimum_dur_sec*1e6;
            %         good_pow_intervals{iC-1} = I(GIX);
            
        end 
    end
    BP = [];


    %         L_gp = L(:,[1 iF+1]);
    ix = (-1*(TS_sfreq/5)):1:TS_sfreq/5;
    PETH_x_sec = ix/TS_sfreq;
    if ~isempty(good_pow_intervals)
        for iN = 1:length(TS)
            if length(TSr{1,iN}) > 20
                %                 figure
                
                %[FR_Burst80,ix,x_sec] = PETH_EEG_simple([bin_ctrs_TS TS_bin(:,iN)],good_pow_intervals{1}(:,1),TS_sfreq/5,TS_sfreq/5,TS_sfreq);
            %    figure
                [FR_Burst80,x_msec] = PETH_raster(TSr{1,iN}/100, good_pow_intervals{1}(:,1)/100,20 , 1000, 1000);

                %                 colorbar
                BP(iN).PETH_M = FR_Burst80;
                BP(iN).PETH_x_sec = x_sec;
                BP(iN).fq_ctrs = fq_ctrs;
            else
                BP(iN).PETH_M = NaN(length(good_pow_intervals{1}), length(PETH_x_sec));
                BP(iN).PETH_x_sec = PETH_x_sec;
                BP(iN).fq_ctrs = fq_ctrs;
            end
        end
    else
        for iN = 1:length(TS)
            BP(iN).PETH_M = nan;
            BP(iN).PETH_x_sec = nan;
            BP(iN).fq_ctrs = fq_ctrs;
        end
    end
 
    
    for iN = 1:length(SP)
        bTBL(tbl_cnt).NeuronID = SP(iN).fname;
        bTBL(tbl_cnt).Session = SES.session;
        bTBL(tbl_cnt).Rat = SES.rat;
        bTBL(tbl_cnt).Interval = ii;
        bTBL(tbl_cnt).Tetrode = SP(iN).Tetrode;
        bTBL(tbl_cnt).Hemisphere = SP(iN).Hemisphere;
        bTBL(tbl_cnt).Depth_uM = SP(iN).Depth_uM;
        bTBL(tbl_cnt).WV = single(SP(iN).WV.mWV);
        
        if length(TSr{iN})> 20
            [a,b]=AutoCorrArray_Abhi(TSr,ac_binsize_ms*1000,n_ac_bins*1000);
        else
            a = nan(1,n_ac_bins);
            b = ac_x_ms;
        end
        bTBL(tbl_cnt).AC = a;
        bTBL(tbl_cnt).AC_x_ms = b/1000;
        
        bTBL(tbl_cnt).fq_ctrs = BP(iN).fq_ctrs;
        bTBL(tbl_cnt).PETH_M = BP(iN).PETH_M;
        bTBL(tbl_cnt).PETH_x_sec = BP(iN).PETH_x_sec;
        bTBL(tbl_cnt).good_pow_intervals = good_pow_intervals;
        
        tbl_cnt = tbl_cnt + 1;
    end
end
    
% TBL = struct2table(pTBL);

OUT.SP = SP;
OUT.bTBL = bTBL;
                
                
%% test plotting stuff
if 0
    for ii = 57:length(pTBL)
        %M = standardize_range(M')';
        %figure
        figure
        subplot(3,1,1:2)
        imagesc(pTBL(ii).PETH_x_sec,[],pTBL(ii).PETH_M)
        %     colormap(jet)
        % Change the color scale.
        caxis_min = prctile(pTBL(ii).PETH_M(~isnan(pTBL(ii).PETH_M)),1);
        caxis_max = prctile(pTBL(ii).PETH_M(~isnan(pTBL(ii).PETH_M)),99);
        caxis([caxis_min caxis_max]);
        pubify_figure_axis
        set(gca,'XTickLabel','')
        plot_vert_line_at_zero
        
        subplot(3,1,3)
        plot_confidence_intervals(pTBL(ii).PETH_x_sec,pTBL(ii).PETH_M)
        % plot(x,nanmean(M))
        axis tight
        box off
        pubify_figure_axis
        plot_vert_line_at_zero
        subplot(3,1,1:2)
    end
end
    
