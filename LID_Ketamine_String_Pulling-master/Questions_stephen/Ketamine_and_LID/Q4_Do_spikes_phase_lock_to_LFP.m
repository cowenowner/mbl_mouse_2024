function OUT = Q4_Do_spikes_phase_lock_to_LFP()
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
% lfp_dir = './LFP';
lfp_dir = DIRS.LFP_Dir;

DO_TREADMILL = true;
n_ac_bins = 100;
ac_binsize_ms = 4;
ac_x_ms = ((1:n_ac_bins)*ac_binsize_ms) - ac_binsize_ms/2;
PLOT_IT = true;
% plot_type = 'wavelet';
SES = LK_Session_Info();
OUT.SES = SES;
OUT.aborted = false;
OUT.n_ac_bins = n_ac_bins;
OUT.ac_binsize_ms = ac_binsize_ms;
tbl_cnt = 1;
OSC = Oscillation_frequencies;
% fq_bands = { 'low_gamma_colgin'  'gamma_80' };

if exist('LFP_INFO.csv','file')
    lfp_info = readtable('LFP_INFO.csv');
else
    lfp_info = readtable('LFP_INFO_auto_generated.csv');
end
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
fq_bands = {'gamma_80'};

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
load('./Processed_Data/best_channels_gamma_80.mat','OUT')
%% Load the best non-reref.
if exist(fullfile('LFP',OUT.best_non_reref),'file')
    % load if stored locally.
    LFP = LK_Load_and_Clean_LFP('./LFP',OUT.best_non_reref);
else
    global DIRS
    lfp_dir = fullfile(DIRS.LFP_Dir,SES.rat_str,SES.session_str,'LFP');
    LFP = LK_Load_and_Clean_LFP(lfp_dir, OUT.best_non_reref);
end

% IMU = LK_Load_and_Process_IMU;
% POS = LK_Load_and_Clean_POS;
% determine what events are available for alignment.
te = []; cnt = 1;
intervals_around_drug_min = [];
if any(E.EventID == 'KetInjectionStart')
    te{cnt} = 'KET';
    e_start_end_uS(cnt,1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
    e_start_end_uS(cnt,2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
    intervals_around_drug_min{cnt} = [-20 -5; 2 30; 30 50];
    
    cnt = cnt + 1;
end
if any(E.EventID == 'LDOPAInjectionStart')
    te{cnt} = 'LDO';
    e_start_end_uS(cnt,1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
    e_start_end_uS(cnt,2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
    intervals_around_drug_min{cnt} = [-20 -5; 2 30; 30 50];
    
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
    
    L = filtfilt(F{1},LFP.LFP(ix1:ix2))';
    t = LFP.t_uS(ix1:ix2);
    [pk,ix] = findpeaks(L);
    
    SP(end).t_uS = t(ix(pk>prctile(pk,60)));
    TS{end} = SP(end).t_uS;
end

pTBL = struct();


for iTE = 1:length(te)
    
    TIMES.(te{iTE}).EventStartEndUsec(1) = e_start_end_uS(iTE,1);
    TIMES.(te{iTE}).EventStartEndUsec(2) = e_start_end_uS(iTE,2);
    TIMES.(te{iTE}).IntervalUsec(1,1) = TIMES.(te{iTE}).EventStartEndUsec(1) + intervals_around_drug_min{iTE}(1,1)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(1,2) = TIMES.(te{iTE}).EventStartEndUsec(1) + intervals_around_drug_min{iTE}(1,2)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(2,1) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(2,1)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(2,2) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(2,2)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(3,1) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(3,1)*60*1e6;
    TIMES.(te{iTE}).IntervalUsec(3,2) = TIMES.(te{iTE}).EventStartEndUsec(2) + intervals_around_drug_min{iTE}(3,2)*60*1e6;
    
    TIMES.(te{iTE}).PeriEventUsec = [intervals_around_drug_min{iTE}(1,1) intervals_around_drug_min{iTE}(2,2) ];
    
    OUT.TM.(te{iTE}) = TIMES;
    
    L = [];
    for iFreq = 1:length(fq_bands)
        for ii = 1:Rows(TIMES.(te{iTE}).IntervalUsec)
            ix1 = binsearch(LFP.t_uS , TIMES.(te{iTE}).IntervalUsec(ii,1));
            ix2 = binsearch(LFP.t_uS , TIMES.(te{iTE}).IntervalUsec(ii,2));
            % Filter data bandpass for each freq and see if spikes phase
            % lock.
            L = filtfilt(F{iFreq},LFP.LFP(ix1:ix2));
            % do something crazy. Look at the response relative to power
            % in the frequency rather than the phase.
            ss = abs(hilbert(L));
            ss = ss - mean(ss);
            Lp = ss;
            
            s = envelope_cowen(abs(L));
            win_sec = 10/mean(OSC.(fq_bands{iFreq})); % Iwant 10 cycles
            win_uS = win_sec*1e6;
            ENV = envelope_cowen(s);
            th = prctile(ENV(1:10:end),90);
            th_l = prctile(ENV(1:10:end),70);
            [SE] = find_intervals([LFP.t_uS(ix1:ix2) ENV], th, th_l, win_uS);
            % create a new subset of spikes.
            % upsample
            TSr = Restrict(TS,SE);
            if isempty(TSr)
                TSr = cell(size(TS));
            end
            
            [PowPhCr]  = SPEC_spike_field_coupling(TSr,[LFP.t_uS(ix1:ix2), Lp]);
            
            if F{iFreq}.HalfPowerFrequency1 >= 40
                % get better phase resolutoin for high frequencies.
                newt_uS = LFP.t_uS(ix1):500:LFP.t_uS(ix2);
                Lu = interp1(LFP.t_uS(ix1:ix2),L,newt_uS,'spline');
                [PhC,phinfo] = SPEC_spike_field_coupling(TS,[newt_uS(:), Lu(:)]);
                [PhCr] = SPEC_spike_field_coupling(TSr,[newt_uS(:), Lu(:)]);
            else
                [PhC,phinfo] = SPEC_spike_field_coupling(TS,[LFP.t_uS(ix1:ix2), L]);
                [PhCr]  = SPEC_spike_field_coupling(TSr,[LFP.t_uS(ix1:ix2), L]);
            end
            
            for iN = 1:length(SP)
                % meta stuff
                pTBL(tbl_cnt).NeuronID = iN;
                pTBL(tbl_cnt).UNeuronID = SES.rat + 1000*SES.session + 100000*iN;
                pTBL(tbl_cnt).Session = SES.session;
                pTBL(tbl_cnt).Rat = SES.rat;
                pTBL(tbl_cnt).Condition = te{iTE};
                pTBL(tbl_cnt).FreqBand = fq_bands{iFreq};
                pTBL(tbl_cnt).Interval = ii;
                pTBL(tbl_cnt).Tetrode = SP(iN).Tetrode;
                pTBL(tbl_cnt).WV = single(SP(iN).WV.mWV);
                pTBL(tbl_cnt).Hemisphere = SP(iN).Hemisphere;
                pTBL(tbl_cnt).Depth_uM = SP(iN).Depth_uM;
                
                if length(TSr{iN})> 20
                    [a,b]=AutoCorr(TSr{iN}/100,ac_binsize_ms,n_ac_bins);
                else
                    a = nan(1,n_ac_bins);
                end
                
                % Calculated stuff
                pTBL(tbl_cnt).Ang = PhC(iN).Ang;
                pTBL(tbl_cnt).hist_rad = PhC(iN).hist_rad;
                pTBL(tbl_cnt).Ang_p = PhC(iN).Ang_p;
                pTBL(tbl_cnt).Ang_to_shuf_p = PhC(iN).Ang_to_shuf_p;
                pTBL(tbl_cnt).Ang_z = PhC(iN).Ang_z;
                pTBL(tbl_cnt).sh_hist_rad_mn = PhC(iN).sh_hist_rad_mn;
                pTBL(tbl_cnt).sh_hist_rad_95ci = PhC(iN).sh_hist_rad_95ci;
                pTBL(tbl_cnt).AC = a;
                
                
                pTBL(tbl_cnt).rAng = PhCr(iN).Ang;
                pTBL(tbl_cnt).rhist_rad = PhCr(iN).hist_rad;
                pTBL(tbl_cnt).rAng_p = PhCr(iN).Ang_p;
                pTBL(tbl_cnt).rAng_to_shuf_p = PhCr(iN).Ang_to_shuf_p;
                pTBL(tbl_cnt).rAng_z = PhCr(iN).Ang_z;
                pTBL(tbl_cnt).rsh_hist_rad_mn = PhCr(iN).sh_hist_rad_mn;
                pTBL(tbl_cnt).rsh_hist_rad_95ci = PhCr(iN).sh_hist_rad_95ci;
                
                
                pTBL(tbl_cnt).powAng = PowPhCr(iN).Ang;
                pTBL(tbl_cnt).powhist_rad = PowPhCr(iN).hist_rad;
                pTBL(tbl_cnt).powAng_p = PowPhCr(iN).Ang_p;
                pTBL(tbl_cnt).powAng_to_shuf_p = PowPhCr(iN).Ang_to_shuf_p;
                pTBL(tbl_cnt).powAng_z = PowPhCr(iN).Ang_z;
                pTBL(tbl_cnt).powsh_hist_rad_mn = PowPhCr(iN).sh_hist_rad_mn;
                pTBL(tbl_cnt).powsh_hist_rad_95ci = PowPhCr(iN).sh_hist_rad_95ci;

                
                tbl_cnt = tbl_cnt + 1;
            end
            
        end
    end
    %
end
TBL = struct2table(pTBL);
TBL.Condition = categorical(TBL.Condition);
TBL.FreqBand = categorical(TBL.FreqBand);
TBL.Hemisphere = categorical(TBL.Hemisphere);
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
        e = envelope_cowen(filtfilt(F{1},LFP.LFP(ix1:ix2)));
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
                        if ii == 1
                            title(['To phase'])
                        end
                    end
                    if ~isempty(v2)
                        subplot(3,3,ii+3)
                        polarhistogram('BinEdges',phinfo.rad_edges,'BinCounts',v2 ,'DisplayStyle','stairs','EdgeColor','k')
                        hold on
                        polarhistogram('BinEdges',phinfo.rad_edges,'BinCounts',TBL.powhist_rad{row_ix},'FaceAlpha',.2,'EdgeColor',clrs(ii,:),'FaceColor',clrs(ii,:))
                        
                        if ii == 1
                            title('To power')
                        end
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

