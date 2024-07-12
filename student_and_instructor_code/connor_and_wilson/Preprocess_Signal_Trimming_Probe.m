%%Reset
clear;
clc;

% Add important paths
addpath(genpath('C:\Users\hp\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\hp\Documents\GitHub\vandermeerlab\code-matlab\photometry'));

% Specify the path of the central folder
centralFolderPath = "C:\Users\hp\Desktop\PhD Stuff\MvdM Rotation\Data\SocialDA\Probe Test";

% Specify the path of the destination folder for saving the results
destinationFolderPath = fullfile(centralFolderPath, 'Trimmed_Signals');
mkdir(destinationFolderPath);

% Get a list of all folders within the central folder
folders = dir(centralFolderPath);
folders = folders([folders.isdir]);  % Select only directories
folders = folders(~ismember({folders.name}, {'.', '..'}));  % Exclude '.' and '..'

% Iterate over each folder
for i = 1:numel(folders)
    folderName = folders(i).name;
    folderPath = fullfile(folders(i).folder, folderName);
    
    % Get the paths of the data files
    cscFilePath = fullfile(folderPath, 'CSC30.ncs');
    eventsFilePath = fullfile(folderPath, 'Events.nev');
    
    

    % Check if the data files exist
    if exist(cscFilePath, 'file') && exist(eventsFilePath, 'file')
        %% Extract and Trim Code
            % -----------------------------------------------------------------------
            
            cd(folderPath)

            % Load experimental keys
      
            LoadExpKeys();
            
            % -----------------------------------------------------------------------
            % manually change paths to these directories 
            addpath(genpath('C:\Users\hp\Documents\GitHub\vandermeerlab\code-matlab\shared'));
            addpath(genpath('C:\Users\hp\Documents\GitHub\vandermeerlab\code-matlab\photometry'));
            
            cfg.fc = {'CSC30.ncs'};
            csc_photo = LoadCSC(cfg);
            
            % extracts FP signal, time, and sampling rate
            FP_data = [];
            FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency; % set FP_data.acq.Fs to sampling frequency rate (5000 points per second) 
            FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1); % set FP_data.acq.Fs time to the time vector subtracted by the first time point
            % this initializes the time vector to start at 0 then 2.0 x 10 ^-4 or 0.0002 and 4.0 x
            %x 10 ^ -4 0.0004 ... 0.0002 second time interval means 5,000 points per
            %second
            FP_data.acq.FP{1} = csc_photo.data'; %fiber data 
            
            % rename variables
            FP = (FP_data.acq.FP{1})*(1239/39); %1239/39 is the voltage multiplier
            time = FP_data.acq.time; 
            FS = FP_data.acq.Fs;
            
            % Start and End Buffer
            % Photobleaching is exponential and often greatest in the first few
            % seconds of recording. This paper recommends removing 2-5 seconds from the
            % beginning and end of the recording file 
            % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7853640/
            
            % Start buffer
            % remove first 2 seconds; number of samples to remove = 2/0.0002;
            FP(1:2/0.0002) = []; 
            time(1:2/0.0002) = [];
            
            % Full Signal Limit
            % We are only interested in the data from Baseline_Start to
            % Conditioning_End. Eveything else can be removed. 
            
            %Slight edits here by Matthew Company
            
            LoadExpKeys();
            cfg_evt = [];
            cfg_evt.eventList = ExpKeys.eventList;
            cfg_evt.eventLabel = ExpKeys.eventLabel;
            evt = LoadEvents(cfg_evt);
            
            % -----------------------------------------------------------------------
            % manually change evt.t events and names 
            
            % Iterate over each cell in the matrix
                for i = 1:numel(evt.t)
                    % Check if the cell contains more than one value
                    if numel(evt.t{i}) > 1
                        % Find the smallest value in the cell
                        minValue = min(evt.t{i});
                        
                        % Keep only the smallest value in the cell
                        evt.t{i} = minValue;
                    end
                end
            


            % Events
            baseline_starts = evt.t{1,2} - csc_photo.tvec(1); % initialze time of events
            baseline_ends = evt.t{1,1} - csc_photo.tvec(1);
            doors_removed = evt.t{1,3} - csc_photo.tvec(1); 
            probe_end = evt.t{1,4} - csc_photo.tvec(1); 
            probe_start = evt.t{1,5} - csc_photo.tvec(1); 

            % Check if evt.t{1,4} exists
            if isfield(evt, 't') && numel(evt.t) >= 4 && ~isempty(evt.t{1,4})
                % evt.t{1,4} exists, do nothing
            else
                % evt.t{1,4} is missing, perform the calculation
                evt.t{1,4} = evt.t{1,5} - 240;
            end
            
            % Check if evt.t{1,3} exists
            if isfield(evt, 't') && numel(evt.t) >= 3 && ~isempty(evt.t{1,3})
                % evt.t{1,3} exists, do nothing
            else
                % evt.t{1,3} is missing, perform the calculation
                evt.t{1,3} = evt.t{1,4} + 1200;
            end


            
            % -----------------------------------------------------------------------
            % manually change variable names based on above. 
            
            % tolerance 
            ind_base_starts = find(abs(time-baseline_starts) < 0.0005);
            ind_base_ends = find(abs(time-baseline_ends) < 0.0005);
            ind_probe_start = find(abs(time-probe_start) < 0.0005);
            ind_doors_removed = find(abs(time-doors_removed) < 0.0005);
            ind_probe_end = find(abs(time-probe_end) < 0.0005);
            
            
            event_time = [time(ind_base_starts), time(ind_base_ends), time(ind_probe_start), time(ind_doors_removed), time(ind_probe_end)];
            event_label = ['Baseline Starts','Baseline Ends','Probe Start','Doors Removed','Probe End'];
            


            %% Downsampled to 1000 Hz 
            % FS is currently 5000 and I'm changing it to 1000 Hz; consider 250 Hz
            dsf = FS/1000;
            FP = decimate(FP,dsf);
            time = downsample(time,dsf);
            FS = FS/dsf;

           %% Raw Signal
            % To plot raw data use FP and multiply by 31.77 (multiplier b/c of voltage
            % divider)
            sessionTitle = 'CW_';
            last_time = length(time); %the value of the last time point is how many seconds the recording was
            timerange1 = 10/0.0002; %datapoint range for 10 s
            timerange2= 100/0.0002; %datapoint range for 100 s 
            timerange3= length(time); % datapoint range for all data
            time_ranges = [timerange1, timerange2, timerange3]; %in seconds 
            
            figure(1)
            for t_i = 1:length(time_ranges)
                t_range = 1:time_ranges(t_i);
                subplot(3, 1, t_i);
                plot(time(t_range), FP(t_range), 'Color', [0 0.5 0])
                title([sessionTitle, num2str(time_ranges(t_i)),'samples'], 'Interpreter','none')
                ylabel('Fiber Signal (V)'); xlabel('Time (s)');
            end
                        
            %% Denoised
            % Why: to filter out electrical noise greater than 10 Hz 
            %       "recording has large electrical noise artifacts, likely due to the
            %       high gain amplifiers in the photodetectors picking up signals from
            %       nearby mobile phone. The artifacts are very short pulses and can be
            %       greatly reduced by running a median filter before the standard low
            %       pass filter. 
            % Method: Median filter & Low pass filter
            % Note: Temporal dynamics of the biosensor are on the level of subseconds
            % (X), so filtering out > 10 Hz signals should be ok. 
            
            % Median filter: remove electrical artifacts 
            FP_denoised = medfilt1(FP);
            
            % check once
            figure(2)
            plot(time,FP);
            title('Median Filtered FP Signal')
            ylabel('Signal (V)')
            xlabel('Time (s)')
            
            % Butterworth Low pass filter 
            % Note: Trisch lab used 20 Hz
            fc = 20; % frequency
            [b,a] = butter(2,fc/(FS/2)); % 2nd order
            %freqz(b,a,[],FP_data.acq.Fs)
            FP_denoised= filter(b,a,FP_denoised);
            xlim([time(1) inf])
            
            figure(3)
            plot(time,FP_denoised)
            title('20 Hz Butterworth Filtered FP Signal')
            ylabel('Signal (V)')
            xlabel('Time (s)')
            xlim([time(1) inf])
            
            figure(4)
            plot(time,FP);
            hold on
            plot(time,FP_denoised)
            hold off
            title('20 Hz Butterworth Filtered FP Signal over Median Filtered FP Signal')
            ylabel('Signal (V)')
            xlabel('Time (s)')
            legend('median','butterworth','Location','northeast')
            xlim([time(1) inf])
            
            %% Detrend 
            % Why: to account for photobleaching 
            % Method: fit an exponential decay to the data and subtract this
            % exponential fit from the signal . 
            % Adam used double exponential fit because of the multiple sources of
            % fluorescence that contributed to the bleaching. (autofluorescence from
            % fiber, brain tissue, and flurophore which may bleach at different rates) 
            
            % This fit is not the best for amphetamine. 
            
            % minimize least square error with lsqcurvefit()
            % input time, FP_denoised 
            t = time;
            y = FP_denoised;
            
            % create a model
            F = @(x,time)x(5)+x(1)*exp(-x(2)*time) + x(3)*exp(-x(4)*time);
            
            % initial parameter guess
            max_sig = max(FP_denoised);
            x0 = [max_sig 0.005 max_sig/2 0.005 0];
            
            % solve least squares
            xunc = lsqcurvefit(F,x0,t,y);
            F_expfit = F(xunc,time);
            
            % subtract the fit from signal 
            FP_detrended = FP_denoised - F_expfit;
            % used for further processing 
            
            figure(5)
            plot(t,FP_denoised)
            hold on
            plot(t,F_expfit)
            hold off
            title('Exponential Fit to Filtered FP Signal')
            ylabel('Signal (V)')
            xlabel('Time (s)')
            xlim([time(1) inf])
            
            figure(6)
            plot(t,FP_detrended)
            title('Detrended and Filtered FP Signal')
            ylabel('Signal (V)')
            xlabel('Time (s)')
            xlim([time(1) inf])

            %% Normalization 
            % Why: to combine data across sessions and/or subjects 
            % "Different sessions may have different levels of fluorophores expression, 
            %  excitation light, and autofluorescence" (Thomas Akam).
            % Method: dF/F or Z-score
            
            % F(t) - F0 / F0; 
            dF = 100.*FP_detrended./F_expfit;
            
            figure(7)
            plot(t,dF)
            title('dF Signal')
            ylabel('Signal dF/F (%)')
            xlabel('Time (s)')
            xlim([time(1) inf])
            
            % Z-score
            % Alternatively, we can normalize by z-scoring each session 
            % subtracting the mean and deviding by standard deviation. 
            % x = current data
            % avg_x = mean of the population ; or alternatively the median 
            % stddev_x = standard deviation of the population 
            
            F_zscored = (FP_detrended - mean(FP_detrended))./std(FP_detrended);
            zdF = (dF - mean(dF))./std(dF);
            
            figure(8)
            plot(t,F_zscored)
            title('Signal z-scored')
            ylabel('Signal (z-scored)')
            xlabel('Time (s)')
            xlim([time(1) inf])
            
            figure(9)
            plot(t,zdF)
            title('Signal dFz-scored')
            ylabel('Signal (dFz-scored)')
            xlabel('Time (s)')
            xlim([time(1) inf])



            %%
            % Initialize variables
            matched_indices = zeros(size(event_time));
            trimmed_F_zscored = [];
            trimmed_t = [];
            trimmed_zdf = [];
            trimmed_df = [];
            
            % Match values from event_time to closest value in t
            for i = 1:numel(event_time)
                [~, idx] = min(abs(t - event_time(i)));
                matched_indices(i) = idx;
            end
            
            % Determine the range for trimming
            trim_start = matched_indices(1,4);
            trim_end = matched_indices(1,5);
            
            % Trim F_zscored and t matrices
            trimmed_F_zscored = F_zscored(trim_start:trim_end);
            trimmed_t = t(trim_start:trim_end);
            trimmed_zdF = zdF(trim_start:trim_end);
            trimmed_dF = dF(trim_start:trim_end);

            
            % Display the matched indices and the trimmed matrices
            disp("Matched Indices:");
            %disp(matched_indices);
            
            disp("Trimmed F_zscored:");
            %disp(trimmed_F_zscored);
            
            disp("Trimmed t:");
            %disp(trimmed_t);

            %%
            figure(10)
            plot(trimmed_t,trimmed_zdF)
            title('Signal dFz-scored')
            ylabel('Signal (dFz-scored)')
            xlabel('Time (s)')
            xlim([trimmed_t(1) inf])


            saveas(gcf, folderName, 'jpg');

            %% Create Save Struct 
            %filename = append(file_name, "processed.mat");
            data.t = trimmed_t;
            data.FP_z = trimmed_F_zscored;
            data.dF = trimmed_dF; 
            data.zdF = trimmed_zdF;
            data.evtt = event_time;
            data.evtlabel = event_label;
            %save(filename, '-struct','data')

            %% Exporting and Saving

    
        % Placeholder code for saving the result as a CSV file
        result = data; %Conversion of final result to result variable
        MATresultFilename = strcat(folderName, '.MAT');
        MATresultFilePath = fullfile(destinationFolderPath, MATresultFilename);
        save(MATresultFilePath, "data")
        
        % End of placeholder code
        
        fprintf('Computation and saving completed for folder: %s\n', folderName);
    else
        fprintf('Missing data files in folder: %s\n', folderName);
    end
end
