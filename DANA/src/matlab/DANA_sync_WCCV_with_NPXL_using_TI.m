function [CV_DATA, WCCV, ALLCV] = DANA_sync_WCCV_with_NPXL_using_TI(WCCV, NPXL_block_start_sec, NPXL_scan_times, TI)
%INPUT NPXL - blocks from neuropixels system
%      WCCV - blcoks from wccv
%      scan_times_NPXL = the timestamps for the start of each scan npxl.
%      OPTIONAL: TI - trial info structure from 
%
% Presume someone fucked up and forgot to start Neuropixels recording
% before FSCV so backtrack from the final scan. NOTE: in the future, just start over.
% Synchronize withe neuropixels data.
PLOT_IT = true;

start_block_npxl = (length(WCCV) - length(NPXL_block_start_sec)) + 1; % sure hope this works.

CV_DATA = [];
st_TI = [TI.start_times_sec];
st_WCCV = [TI.start_times_sec];
Hz_TI = round([TI.Hz]);
LV_TI = [TI.LV];
HZ_WCCV = [WCCV.Hz];
LV_WCCV = [WCCV.LV];
% Pad in case not aligned.

if PLOT_IT
figure
subplot(1,3,1)
plot(HZ_WCCV,'o-')
hold on
plot(Hz_TI,'*r-')
ylabel('Hz')

xcorr(Hz_TI,HZ_WCCV,10)


subplot(1,3,2)
plot(LV_WCCV,'o-')
hold on
plot(LV_TI,'*r-')
ylabel('LV')

subplot(1,3,3)
plot(st_TI,st_WCCV,'o')
yyaxis right
plot(pdist([st_TI st_WCCV]))

figure
[CC,lags] = xcorr(LV_TI,LV_WCCV,11)
plot(lags,CC);plot_vert_line_at_zero;
[~,ix] = max(CC);
shft = lags(ix);
end

for iBlock = 1:length(WCCV)
    % Find the closest match in TI
    % [min_val, closest_time_ix] = min(abs(st_WCCV(iBlock)-st_TI));
    % if min_val > .1
    %     continue
    % end
    t = Restrict(NPXL_scan_times, NPXL_block_start_sec(iBlock)-.001, NPXL_block_start_sec(iBlock) + WCCV(iBlock).block_dur_sec + .001);
    fprintf('%d) WCCV %d scans, NPXL %d scans\n',iBlock,Rows(WCCV(iBlock).within_block_CV_data),length(t));
    % assign the data to these points.
    CV_DATA = [CV_DATA; t(:) WCCV(iBlock).within_block_CV_data(1:length(t),2)];
    WCCV(iBlock).NPXL_timestamps_sec = t(:);
end
if nargout > 2
    % Storing the data in this format can often be more useful.
    ALLCV = [];
    ALLCV.M = [];
    ALLCV.Hz = [];
    ALLCV.LV = [];
    ALLCV.x_sec = WCCV(1).within_block_CV_data(:,1)';

    for iB = 1:length(WCCV)
        ALLCV.M(iB,:) = WCCV(iB).within_block_CV_data(:,2)';
        ALLCV.Hz(iB) = WCCV(iB).Hz;
        ALLCV.LV(iB) = WCCV(iB).LV;
        ALLCV.Start_sec(iB) = WCCV(iB).NPXL_timestamps_sec(1);
    end
end

if nargin == 0
    figure
    plot(CV_DATA(:,1), CV_DATA(:,2))
end
