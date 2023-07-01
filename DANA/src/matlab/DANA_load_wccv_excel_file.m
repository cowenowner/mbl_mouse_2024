function [DATA, STIM] = DANA_load_wccv_excel_file(WCCV_path_and_fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the CV data and stimulation times from the excel sheet.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%
upsample_sample_rate_Hz = 20;
CV_smooth_win_sec = .5;
stim_fq_string = {'10Hz_vs_time' '20Hz_vs_time'};
output_vbl = 'Concentration';
win_bins = round(CV_smooth_win_sec*upsample_sample_rate_Hz);
han_win = hanning(win_bins)/sum(hanning(win_bins));
time_to_add_to_stim_start_sec = 0; % this is in case there is a systematic difference (Which there should not be) between the stim times data and the dopamine measurement times.
DATA = [];
STIM = [];
[~,fname] = fileparts(WCCV_path_and_fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the stimulation times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = readtable(WCCV_path_and_fname,'Sheet', 'Stim_trains','NumHeaderLines',1);

for iC = 1:length(T.Properties.VariableNames)
    [str,rem] = strtok(T.Properties.VariableNames{iC},'_');
    if ~isempty(rem)
        STIM(iC).Hz = str2double(str(2:end-2));
        STIM(iC).LV = str2double(rem(4)) + str2double(rem(end-1:end))/100 ;
        STIM(iC).Stim_time_s = T.(T.Properties.VariableNames{iC}) + time_to_add_to_stim_start_sec;
        [STIM(iC).Stim_time_IFR, STIM(iC).Stim_time_IFRsmth, t] = ...
            Instantaneous_Firing_Rate(STIM(iC).Stim_time_s,0.02,40,'start_time',0,'PLOT_IT',false);
        STIM(iC).Stim_time_IFR_t_sec = t(1:end-1) + 0.02/2; % center the bins
    end
end
% Load the CV data...
% Probably best to create a new table for each frequency That is how the
% data is stored. Keeping close to the data might make it easier...
cnt = 1;
for iHz = 1:length(stim_fq_string)
    fq = str2double(stim_fq_string{iHz}(1:2));
    T = readtable(WCCV_path_and_fname,'Sheet', stim_fq_string{iHz},'NumHeaderLines',2 );
    if isempty(T)
        continue
    end
    if isnan(T.Time_s(end))
        error('There is a missing timestamp (last one) in the excel file. Fix this.')
    end

    LV_row = readtable(WCCV_path_and_fname,'Sheet', stim_fq_string{iHz},'Range', 'A2:CF2');
    upsample_T_sec = T.Time_s(1):1/upsample_sample_rate_Hz:T.Time_s(end);
    upsample_T_sec = upsample_T_sec(:);
    % Go through each column and store as an element in a structure array
    for iCol = 2:length(T.Properties.VariableNames)
        str = T.Properties.VariableNames{iCol};
        if contains(str,output_vbl)

            LV = table2array(LV_row(1,iCol));
            [v,rem] = strtok(str,'_');
            if ~isempty(rem)
                ix = find([STIM.Hz] == fq & [STIM.LV] == LV);
                if isempty(ix) || isempty(STIM(ix).Stim_time_s)
                    disp('Error in reading in the excel file')
                    continue
                end

                DATA(cnt).Stimulation_Hz = fq;
                DATA(cnt).LV = LV;
                DATA(cnt).Time_s = T.Time_s;
                DATA(cnt).Data_type = v;
                DATA(cnt).Data = T.(str);
                DATA(cnt).Trial = str2double(rem(2:end));

                DATA(cnt).upsample_Time_s = upsample_T_sec;
                v = interp1(DATA(cnt).Time_s,DATA(cnt).Data, upsample_T_sec);
                v = conv(v,han_win,'same');
                DATA(cnt).upsample_Data = v(:);

                DATA(cnt).Stim_time_s = STIM(ix).Stim_time_s;
                DATA(cnt).Stim_time_IFRsmth = STIM(ix).Stim_time_IFRsmth;
                DATA(cnt).Stim_time_IFR_t_sec = STIM(ix).Stim_time_IFR_t_sec;
                [DATA(cnt).Stim_sliding_LV,DATA(cnt).Stim_sliding_LVR, INFO] = LocalVariance(diff(STIM(ix).Stim_time_s*1000),[],[4000 1000]);
                DATA(cnt).Stim_sliding_LV_t_sec = mean(INFO.start_and_end_times/1000,2)+STIM(ix).Stim_time_s(1);
                cnt = cnt + 1;
            end
        end
    end
end


if nargout == 0
    %% Plot
    LV = 1; Hz = 10;
    %     LV = 1.24;
    IX =  [DATA.LV] == LV & contains({DATA.Data_type},'Concentration') & [DATA.Stimulation_Hz] == Hz;
    ix = find(IX);
    M = [DATA(IX).Data];
    TIME_IX = DATA(ix(1)).Time_s < 30;
    ss = DATA(ix(1)).Stim_time_s(1);
    figure
    subplot(6,1,1:2)
    imagesc(DATA(ix(1)).Time_s(TIME_IX)-ss,[],M(TIME_IX,:)');
    colorbar_label('[DA]nM')
    hold on
    plot_tics(DATA(ix(1)).Stim_time_s-ss,[1 0 0],1,1)
    title(sprintf('%s LV %1.2f, %1.1f Hz Concentration',fname,LV,Hz))
    plot_vert_line_at_zero;

    subplot(6,1,3:4)
    plot_confidence_intervals(DATA(ix(1)).Time_s(TIME_IX)-ss,M(TIME_IX,:)')
    ylabel('[DA]')

    hold on
    plot_tics(DATA(ix(1)).Stim_time_s-ss,[1 0 0],0,50)
    yyaxis right
    plot(DATA(ix(1)).Stim_time_IFR_t_sec-ss,DATA(ix(1)).Stim_time_IFRsmth,'b')
    xlabel('sec')
    ylabel('[DA]')
    plot_vert_line_at_zero;
    ylabel('IFR Hz')

    subplot(6,1,5:6)
    Md = diff(M);
    Md = [Md;Md(end,:)];
    plot_confidence_intervals(DATA(ix(1)).Time_s(TIME_IX)-ss,Md(TIME_IX,:)')
    hold on
    plot(DATA(ix(1)).Stim_time_s-ss,zeros(size(DATA(ix(1)).Stim_time_s)),'r+')
    ylabel('delta [DA]')
    yyaxis right
    plot(DATA(ix(1)).Stim_time_IFR_t_sec-ss,DATA(ix(1)).Stim_time_IFRsmth,'b')
    xlabel('sec')
    plot_vert_line_at_zero;

    %%%%%%%%%%%%%%
    % Inter-trial correlations
    STIM_TIME_IX = DATA(ix(1)).Time_s > 10 & DATA(ix(1)).Time_s < 20;
    figure
    subplot(2,2,1)
    imagesc(corr(M(STIM_TIME_IX,:)))
    axis square
    xlabel('Trial'); ylabel('Trial')
    caxis([.01 .99])
    colorbar
    title('[DA] Pearson r just during stim')
    subplot(2,2,2)
    imagesc(corr(Md(STIM_TIME_IX,:)))
    axis square
    caxis([.01 .99])
    colorbar
    title('delta[DA]')
    % Noise correlations...
    % Inter-trial correlations

    subplot(2,2,3)
    imagesc(corr(M(STIM_TIME_IX,:) - mean(M(STIM_TIME_IX,:),2)))
    axis square
    xlabel('Trial'); ylabel('Trial')
    caxis([-.5 .5])
    colorbar
    title('noise corr Pearson r')
    subplot(2,2,4)
    imagesc(corr(Md(STIM_TIME_IX,:) - mean(Md(STIM_TIME_IX,:),2)))
    axis square
    caxis([-.5 .5])
    colorbar

    % Scatter plots
    figure

    for iC = 1:Cols(M)
        Mup(:,iC) = interp1(DATA(ix(1)).Time_s(STIM_TIME_IX), M(STIM_TIME_IX,iC),DATA(ix(1)).Stim_time_IFRsmth);
        Mup_delta(:,iC) = interp1(DATA(ix(1)).Time_s(STIM_TIME_IX), Md(STIM_TIME_IX,iC),DATA(ix(1)).Stim_time_IFRsmth);
    end
    subplot(2,2,1)
    for iC = 1:Cols(M)
        scatter(DATA(ix(1)).Stim_time_IFRsmth,Mup(:,iC))
        hold on
        lsline
    end
    xlabel('IFR');ylabel('[DA]')
    subplot(2,2,2)
    for iC = 1:Cols(Md)
        scatter(DATA(ix(1)).Stim_time_IFRsmth,Mup_delta(:,iC))
        hold on
        lsline
    end
    xlabel('IFR');ylabel('d[DA]')

    subplot(2,2,3)
    for iC = 1:Cols(Md)
        D = [DATA(ix(1)).Stim_time_IFRsmth Mup(:,iC)];
        GIX = ~isnan(sum(D,2));
        D = D(GIX,:);
        [c,lags] = xcorr(D(:,1),D(:,2),'coeff');
        plot(lags,c)
        hold on
    end
    axis tight
    xlabel('lag (s)')
    ylabel('r')
    title('xcorr stim to DA')
    plot_vert_line_at_zero

    subplot(2,2,4)
    for iC = 1:Cols(Md)
        D = [DATA(ix(1)).Stim_time_IFRsmth Mup_delta(:,iC)];
        GIX = ~isnan(sum(D,2));
        D = D(GIX,:);
        [c,lags] = xcorr(D(:,1),D(:,2),'coeff');
        plot(lags,c)
        hold on
    end
    axis tight
    xlabel('lag (s)')
    ylabel('r')
    plot_vert_line_at_zero
    title('xcorr stim to delta DA')

    %%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    TIME_IX = DATA(ix(1)).upsample_Time_s < 30;
    M = [DATA(IX).upsample_Data];
    subplot(4,1,1:2)
    imagesc(DATA(ix(1)).upsample_Time_s(TIME_IX)-ss,[],M(TIME_IX,:)');
    hold on
    plot_tics(DATA(ix(1)).Stim_time_s-ss,[1 0 0],1,1)
    plot_vert_line_at_zero;

    title(sprintf('%s LV %1.2f, %1.1f Hz Concentration',fname,LV,Hz))
    subplot(4,1,3:4)
    plot_confidence_intervals(DATA(ix(1)).upsample_Time_s(TIME_IX)-ss,M(TIME_IX,:)')
    hold on
    % plot(DATA(ix(1)).Stim_time_s-ss,zeros(size(DATA(ix(1)).Stim_time_s)),'r+')
    plot_tics(DATA(ix(1)).Stim_time_s-ss,[1 0 0],0,50)
    xlabel('sec')
    ylabel('[DA]')
    plot_vert_line_at_zero;

    % STAs
    M = [DATA(IX).Data];
    PETH_spike_triggered_average([DATA(ix(1)).Time_s*1000 M(:,1)], DATA(ix(1)).Stim_time_s*1000, 2000, 2000, 5, 1000, true)

end