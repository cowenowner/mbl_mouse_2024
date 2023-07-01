function [stim_train,CV_train,blank_intervals]= DANA_eliminate_overlap_with_CV_pulses(stim_train,CV_freq, CV_duration_s, stim_pulse_dur_sec, min_isi_s, stim_train_dur_s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure that there is no overlap between the CV pulses and the
% individual stim pulses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The CV pulse train.
% TODO: Ensure that the first spike is just after the first CV pulse and
% that the last is just before the last CV pulse. That way all trains have
% the exact same duration. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    CV_freq = 5; % once every 200 ms.
end
if nargin < 3
    CV_duration_s = 0.0088; % the typical duration of the ramp up/down of the CV pulse.
end
if nargin < 4
    stim_pulse_dur_sec = 0.005; % the duration of the stimulation pulse so that it does not overlap with the CV interval
end
if nargin < 5
    min_isi_s = 0.005; % minimum time between stimulations
end
if nargin < 6
    stim_train_dur_s = stim_train(end);
end
% Get rid of stim pulses that are just too close together.
% Maybe not - how about just increase the ISI for that pulse and subtract
% it from the subsequent ISI.
stim_train = stim_train(:);
orig_stim_train = stim_train;

ISIs_orig = diff(stim_train);
ISIs = ISIs_orig;
BIX = ISIs_orig<min_isi_s;
ISIs(BIX) = ISIs_orig(BIX) + min_isi_s;
if any(BIX)
    sprintf('Found %d small inter-stim intervals. Removed', sum(BIX))
end
% To preserve the mean firing rate, subtract the added amount of time from
% the overall ISI array.
added_time = sum(ISIs) - sum(ISIs_orig);
if added_time > 0
    bigger_ISIs_IX = ISIs>median(ISIs);
    ISIs(bigger_ISIs_IX) = ISIs(bigger_ISIs_IX) - added_time/sum(bigger_ISIs_IX);
end
% Just in case the above made a few ISIs get too small...
BIX = ISIs<min_isi_s;
ISIs(BIX) = ISIs(BIX) + min_isi_s;


% update the stim_train
stim_train = [stim_train(1); stim_train(1) + cumsum(ISIs)];

CV_train = 0:1/CV_freq:max(stim_train);
CV_train = CV_train(:);
n_shifts = 0;
if ~isempty(CV_duration_s)
    % intervals for elimination of pulses as they overlap with CV pulses.
    % subtract stim_params.stim_pulse_dur_sec so that a pulse does not
    % bleed into the CV pulse.
    blank_intervals = [CV_train(:)-stim_pulse_dur_sec CV_train(:) + CV_duration_s];
    for ii = 1:Rows(blank_intervals)
        IX = stim_train >= blank_intervals(ii,1) & stim_train <= blank_intervals(ii,2);
        if any(IX)
            n_shifts = n_shifts + sum(IX);
            % just move them to the end of the CV pulse to that they do not overlap.
            % seems reasonable since this is a rare event.
            % NOTE: this is only one direction so the mean effect would be
            % to lower the frequency of firing. A smarter way would be to
            % also shorten the previous ISI
            % Could just randomize - sometimes add, sometimes subtract.
            % That would keep the order preserved - could have some
            % timestamps go before others though
            % Assign it to the nearest interval.
%             ix = find(IX);
            % this could create some new problems for high frequency stim-
%             % will get multiple stims - should still be rare though.
%             for jj = 1:length(ix)
%                 stim_train(ix(jj)) = Closest(blank_intervals(:),stim_train(ix(jj)));
%             end
            % the following seems reasonable as well.
             if rand(1,1)>.5
                 stim_train(IX) = stim_train(IX) + CV_duration_s + stim_pulse_dur_sec;
             else
                 stim_train(IX) = stim_train(IX) - CV_duration_s - stim_pulse_dur_sec;
             end
        end
    end
    disp(['Shifted ' num2str(n_shifts) ' stims']);
end
stim_train = unique(stim_train); % sorts them in the case that things went out of order.

% The above manipulations might have caused some ISI issues
ISIs = diff(stim_train);
BIX1 = ISIs<min_isi_s;
BIX2 = true;
% jitter stuff around until there is no overlap
n_loops = 1;
stim_train1 = stim_train;
while any(BIX1) || any(BIX2)
    ISIs(BIX1) = min_isi_s + rand(sum(BIX1),1) * 4 * min_isi_s;
    stim_train = [stim_train(1); stim_train(1) + cumsum(ISIs)];
    stim_train = sort(stim_train);
    stim_train(end) = stim_train_dur_s; % enforce that there is a spike at the start and end of the sequence.
    % Blank intervals
    [~,BIX2] = Restrict(stim_train,blank_intervals);
    stim_train(BIX2) = stim_train(BIX2) + .058*rand(sum(BIX2),1) - .058/2;
    % before 0
    stim_train(stim_train<0)= blank_intervals(1,2) + .1*rand(sum(stim_train<0),1);
    stim_train = sort(stim_train);

    % after last time.
    IX = (stim_train>stim_train_dur_s);
    stim_train(IX) = blank_intervals(end-2,2) + .1* rand(sum(IX),1);
    stim_train = sort(stim_train);
    % in blank intervals
    [~,BIX2] = Restrict(stim_train,blank_intervals);
    stim_train(BIX2) = stim_train(BIX2) + .5*rand(sum(BIX2),1)-.25; % This seems WAY to big. I think
    stim_train = sort(stim_train);
    
    tmpISI = diff(stim_train);
    tmpISI(tmpISI<min_isi_s) = min_isi_s + .05*rand(1,1);
    stim_train = [stim_train(1); stim_train(1) + cumsum(tmpISI)];

    BIX1 = diff(stim_train)<min_isi_s;
    [~,BIX22] = Restrict(stim_train,blank_intervals);
    BIX2 = stim_train<0 | (stim_train>stim_train_dur_s) | BIX22;

    n_loops = n_loops + 1;
     if n_loops > 20000
         disp('wtf')
%          stim_train = stim_train1 + randn(size(stim_train1))*.001;
%          break
     end
end
n_loops
% double check -should never need this but to be safe...
if ~isempty(CV_duration_s)
    for ii = 1:Rows(blank_intervals)
        IX = stim_train >= blank_intervals(ii,1) & stim_train <= blank_intervals(ii,2);
        if any(IX)
            error('this should never happen. found overlap')
        end
    end
end
