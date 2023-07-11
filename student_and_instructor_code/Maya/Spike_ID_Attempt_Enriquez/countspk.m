function [totalspk] = countspk(start_t, end_t, SP) % Define output/inputs 
    SP.t_sec = (0:(Cols(SP)-1))/1e-6; %Convert microseconds to seconds
    totalspk = zeros(size(SP)); % Make empty variable same length of dataset

    for i = 1:length(SP) %For every timepoint in the spike time dataset
     output =  SP(i).t_sec > start_t & SP(i).t_sec < end_t; %find the spikes between these timepoints
     totalspk(i) = sum(output); % Add up the spikes within these timepoints
    end
totalspk = totalspk'; %Invert so is column, not row   

end