function E = Events_from_clock_and_signal_times(clock_times_sec, signal_times_sec, min_interval_between_packets_sec, code_length)
% Assumes packets (code_length) are 4 bits long unless specified otherwise.
%
% 
% Cowen 2023
if nargin < 3
    min_interval_between_packets_sec = min(diff(clock_times_sec))*1.3;
end
if nargin < 4
    code_length = 4; % number of bits in the code
end

%%
interval_between_clock_sec = min(diff(clock_times_sec));

diff_sec = diff([0; clock_times_sec(:)]);


code_ID = zeros(size(clock_times_sec));
dist_sec = zeros(size(clock_times_sec));
code_time_sec = zeros(size(clock_times_sec));

code_cnt = 1;
event_time_sec = clock_times_sec(1);
for ii = 1:length(clock_times_sec)

    ix = Closest(signal_times_sec,clock_times_sec(ii));
    dist_sec(ii) = abs(clock_times_sec(ii) - signal_times_sec(ix));

    if diff_sec(ii) > min_interval_between_packets_sec
        code_cnt = code_cnt + 1;
        event_time_sec = clock_times_sec(ii);
    end
    code_ID(ii) = code_cnt;
    code_time_sec(ii) = event_time_sec; % time of the start of the code.

end
code_val = dist_sec < interval_between_clock_sec/2;

uCode = unique(code_ID);
code_cnt = 1;
for iU = 1:length(uCode)
    ix = find(code_ID == uCode(iU));
    BIN = code_val(ix)';
    
    if length(BIN) ~= code_length

        if length(BIN) == code_length*2
            % we have 2 packets merged into one. Split them.
            E.Event_time_sec(code_cnt) = code_time_sec(ix(1));
            E.Event_code(code_cnt) = bin2dec_matrix(BIN(1:(code_length)));
            code_cnt = code_cnt + 1;
            E.Event_time_sec(code_cnt) = code_time_sec(ix(1)) + .001;
            E.Event_code(code_cnt) = bin2dec_matrix(BIN((code_length*2-1):(code_length*2)));
            code_cnt = code_cnt + 1;
        elseif length(BIN) == code_length*3
            % we have 3 packets merged into one. Split them.
            E.Event_time_sec(code_cnt) = code_time_sec(ix(1));
            E.Event_code(code_cnt) = bin2dec_matrix(BIN(1:(code_length)));
            code_cnt = code_cnt + 1;
            E.Event_time_sec(code_cnt) = code_time_sec(ix(1)) + .001;
            E.Event_code(code_cnt) = bin2dec_matrix(BIN((code_length*2-1):(code_length*2)));
            code_cnt = code_cnt + 1;
            E.Event_time_sec(code_cnt) = code_time_sec(ix(1)) + .001;
            E.Event_code(code_cnt) = bin2dec_matrix(BIN((code_length*3-1):(code_length*3)));
            code_cnt = code_cnt + 1;

        else
            % some strange packet.
            fprintf('BAD: code exeeds length code %d len %d \n',iU, length(BIN))

            E.Event_time_sec(code_cnt) = nan;
            E.Event_code(code_cnt) = nan;
            code_cnt = code_cnt + 1;

        end
    else

        E.Event_time_sec(code_cnt) = code_time_sec(ix(1));
        E.Event_code(code_cnt) = bin2dec_matrix(BIN);      
        code_cnt = code_cnt + 1;
    end
end
if nargout == 0
     
    figure
    plot(clock_times_sec,zeros(size(clock_times_sec)),'.')
    hold on
    plot(signal_times_sec,zeros(size(signal_times_sec)),'ro')

end
