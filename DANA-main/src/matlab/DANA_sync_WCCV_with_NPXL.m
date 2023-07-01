function [CV_DATA, WCCV, ALLCV] = DANA_sync_WCCV_with_NPXL(WCCV, NPXL_block_start_sec, NPXL_scan_times, TI)
%INPUT NPXL - blocks from neuropixels system
%      WCCV - blcoks from wccv
%      scan_times_NPXL = the timestamps for the start of each scan npxl.
%      OPTIONAL: TI - trial info structure from 
%
% Presume someone fucked up and forgot to start Neuropixels recording
% before FSCV so backtrack from the final scan. NOTE: in the future, just start over.
% Synchronize withe neuropixels data.
if nargin < 4
    TI = [];
end
start_block_npxl = (length(WCCV) - length(NPXL_block_start_sec)) + 1; % sure hope this works.
if start_block_npxl == -1
    start_block_npxl = 1;
end
WCCV = WCCV(start_block_npxl:end); % assume WCCV was recording before NPXL started so reduce it.
if length(WCCV) ~= length( NPXL_block_start_sec )
    error('mismatch block times.')
end
CV_DATA = [];
for iBlock = 1:length(WCCV)
    % Find the start and end of this block in neuropixels times.
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
    plot(CV_DATA(:,1), CV_DATA(:,3))
end
