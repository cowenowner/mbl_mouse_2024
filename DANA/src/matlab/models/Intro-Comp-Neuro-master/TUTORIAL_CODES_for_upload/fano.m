function [ fano_factor, time_windows ] = fano( spikes, dt, method, window_size )
%fano( spikes, dt, method, window_size ) produces the Fano factor for a series of events
%   spikes is an array of 1s and 0s with a 1 indicating the event
%   If spikes is a vector (1D array) the function will split the vector
%   into different size time bins and calculate the variance divided by the
%   mean of the number of events in each time bin.
%   If spikes is a two-dimensional array, the successive rows contain
%   successive trials. The Fano factor can be evaluated by counting the
%   number of spikes in the time interval from the start of each trial to
%   each time point then comparing those cumulative sums across trials 
%   (method 1) or by taking a fixed size of time-window and shifting it 
%   across the duration of the trial (method 2).
%
%   dt is the time interval between successive columns in the spikes array.
%
%   method is used if ndims(spikes) == 2.
%   method can be 1 or 2 (see below). Default method=1.
%
%   window_size can be set and is used as a fixed value for method 2.
%   By default window_size = 100ms.
%
%   If spikes are from a single time series they are in an array of rank 1.
%   In this case we split up the time series into bins of increasing size
%   and take the average number and variance in number of spikes across all
%   bins in the one time series.
%
%   If spikes are from many trials, then in method 1, the time window
%   commences at the first time point of the spike array and extends
%   through the complete trial, so spikes up to a given point are
%   accumulated using "cumsum". Alternatively, in method 2, a default
%   window size is used (that can be chosen at input) and the time-window
%   slides across the trial.
%
%   This function is required for Tutorial 3.2 in the textbook
%   An Introductory Course in Computational Neuroscience
%   by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('method','var') )  % If not sent to the function
    method = 1;         % use the default value of 1
end
if (~exist('window_size','var') ) % If not sent to the function
    window_size = 0.1;      % use the default value of 0.1
end
nwindow_size = round(window_size/dt);   % Number of timebins in fixed window

if ( ndims(spikes) > 2)  
            disp('ERROR! too many dimensions in spikes array')
        fano_factor = [];
        time_windows = [];
else
    if ( min(size(spikes)) == 1 )               % If spikes is a 1D array    
        Nt = length(spikes);                    % Number of time_bins
        window_min = round(0.010/dt);
        time_windows = window_min:window_min:round(Nt/25);       % Different sizes to loop over.
        fano_factor = zeros(size(time_windows));   % Initialize to zero
        ff_bin = 1;                                % Index for the array
        for nwindow_size = time_windows;           % Now loop through bin-sizes
            num_windows = floor(Nt/nwindow_size);  % Number of windows of this size
            N = zeros(1,num_windows);              % No. of sikes per bin
            for i = 1:num_windows;                 % For each window, count spikes
                N(i) = sum(spikes((i-1)*nwindow_size+1:i*nwindow_size));
            end
            fano_factor(ff_bin) = var(N)/mean(N);  % Calculate Fano factor for that bin-size
            ff_bin = ff_bin + 1;                   % Increment index for next calculation
        end
    else                            % If there are many trials so spikes is 2D
        if ( method == 1 )                      % Default, cumulative trial method
            Nspikes = cumsum(spikes,2);         % Cumulative sum of spikes in each trial
            Nvar = var(Nspikes,0,1);            % Variance across trials
            Nmean = mean(Nspikes,1);            % Mean across trials
            fano_factor = Nvar./Nmean;          % Fano factor is variance/mean
            time_windows = 1:length(fano_factor); % Indices for time_windows
        else
            if ( method == 2 )                  % Alternative method with fixed width sliding window
                [~, Nt] = size(spikes);         % Obtain number of time bins 
                Nt = Nt + 1 - nwindow_size;     % Alter to be number of windows
                time_windows = round(nwindow_size/2):Nt-1+round(nwindow_size/2);
                fano_factor = zeros(1,Nt);      % Initialize output vector
                for i = 1:Nt;                   % Loop through time bins
                    % Sum across time index to find no. of spikes in each bin 
                    N = sum(spikes(:,i:i+nwindow_size-1),2);   
                    fano_factor(i) = var(N)/mean(N);    % Calculate Fano factor
                end
            else
                disp('ERROR! method must be 1 or 2')
                fano_factor = [];
                time_windows = [];
            end
        end
    end
end

time_windows = dt*time_windows;     % Return time_windows in units of time

end

