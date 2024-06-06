function [fixed_sync_times] = Fix_missing_video_sync_pulses(sync_times, bad_interval_threshold)
% Sometimes video sync pulses can go missing - for example if the cable
% gets shorted or disconnected briefly. This code is inteded to help fill
% in relatively short intervals where TTL pulse data is missing.
%
% In some cases, the video will be turned on and off in blocks. The user
% will be asked to identify the blocks Then, within each block, the code
% will look for missing pulses and if found, will fill in the missing
% pulses with its best guess.
%
% The user will be prompted to select a time ranges that constitued
% contiguous good blocks of events. Then, within each block, missing
% timestmaps will be 'filled in'.
%
% INPUT:
% list of timestamps (you choose the units)
% the expected inter-pulse interval. Anything 1.5x bigger than this in the
% data will be considered data elegible for 'filling in'.
%
% OUTPUT
% cleaned up timestmamps - by filling in missing data with data using the
% same inter-pulse interval
%
% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix_missing_video_sync_pulses(EVT.frontcam.time_sec, 0.00283200000012584)
if nargin < 2
    bad_interval_threshold = median(diff(sync_times))*1.5;
end

bad_ix = find(diff(sync_times)>bad_interval_threshold);
df = diff(sync_times);
md_diff = median(df);
fixed_sync_times = [];
figure
subplot(2,2,1:2)
plot(sync_times,zeros(size(sync_times)),'.','MarkerSize',1)
hold on
plot(sync_times(bad_ix),zeros(size(bad_ix)),'ro')
title(sprintf('median inter interval %d, min %d, max %d',md_diff, min(df), max(df)))

subplot(2,2,3)
histogram(diff(sync_times))
xlabel('inter pulse intervals')
subplot(2,2,4)
histogram(log10(diff(sync_times)))
xlabel('log inter pulse intervals')

answer = questdlg('Do the timestamps look correct?','Answer', 'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
        all_OK = true;
    case 'No'
        all_OK = false;
end
if all_OK
    % everything is fine, just return the original timestamps.
    fixed_sync_times = sync_times;
else
    % Something ain't right
    answer = inputdlg({'How many good blocks are there?'},'Input',[1 35], {'1'});
    n_blocks = str2double(answer{1});

    figure(1)
    clf
    plot(sync_times,zeros(size(sync_times)),'.','MarkerSize',5)
    hold on
    plot(sync_times(bad_ix),zeros(size(bad_ix)),'ro')
    new_times = [];
    for iInterval = 1:n_blocks
        title('Click at the start and end of each good block (do not miss the start or end)')

        x = ginput(2);
        x = x(:,1);
        plot(x(1),0,'g>')
        plot(x(2),0,'r<')
        %%%%%%%%%%
        new_times = sync_times(sync_times >= x(1) & sync_times <= x(2));
        diffs = diff(new_times(:))';
        bad_ix2 = find(diffs>bad_interval_threshold);
        % for each of these intervals, fill in with good data.
        st_ed_ix = [bad_ix2(:) bad_ix2(:)+1];
        st_ed_time = new_times(st_ed_ix);
        fill_times = [];
        md_diff = median(diffs(diffs<bad_interval_threshold));
        for iR = 1:Rows(st_ed_time)
            % fill_times = [fill_times (st_ed_time(iR,1)+md_diff):md_diff:(st_ed_time(iR,2)-md_diff)];
            se = [st_ed_time(iR,1)+md_diff st_ed_time(iR,2)-md_diff];
            dur = diff(se);
            fill_times = [fill_times linspace(se(1), se(2),ceil(dur/md_diff)+1)];
            
            if st_ed_time(iR,2) - fill_times(end) > bad_interval_threshold *.95;
                fill_times(end+1) = fill_times(end) + md_diff*.95;
            end
        end
        updated_times = unique([fill_times(:); new_times(:)]);
        big_ix = find(diffs > bad_interval_threshold);
        small_ix = find(diffs < md_diff/2.5);
        updated_times(small_ix) = [];
        if 0
            figure
            plot(new_times,zeros(size(new_times)),'+k')
            hold on
            plot(updated_times,zeros(size(updated_times)),'or')
            plot_markers_simple(st_ed_time)
            plot(updated_times(big_ix), zeros(size(big_ix)),'c*')
            min(diff(updated_times))
            max(diff(updated_times))
            diffs = diff(updated_times);
        end
        fixed_sync_times = [fixed_sync_times; updated_times(:)];
        % plot(st_ed_time(:,1),zeros(size(st_ed_time(:,1)))+.1,'>g')
        % plot(st_ed_time(:,2),zeros(size(st_ed_time(:,2)))+.1,'<r')
        % st_ed_time(:,2) - st_ed_time(:,1);
        % find intervals with missing data.
        % figure(2)
        % histogram(diffs(diffs<= median(diffs)*1.1))
    end

end
