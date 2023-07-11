%% Identifying Spikes
% This script is the first attempt to make a function that will count the
% spikes in each channel within a specific timepoint of spike dataset SP

%% Creating variable for timepoint 

%Test with a single channel
%output =  SP(1).t_uS > 0.5e+9 & SP(1).t_uS < 1e+9;
%totalspk = sum(output)

% Loop for finding all the spikes within a certain timepoint

totalspk = zeros(size(SP)); %Make empty variable same length of dataset
for i = 1:length(totalspk) %For every timepoint in the spike time dataset
    output =  SP(i).t_uS > 0.5e+9 & SP(i).t_uS < 1e+9; %find the spikes between these timepoints
    totalspk(i) = sum(output); % Add up the spikes within these timepoints
end
totalspk = totalspk'; %Invert so is column, not row 



