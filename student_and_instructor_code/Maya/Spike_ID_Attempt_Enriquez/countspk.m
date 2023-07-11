function [totalspk] = countspk(start_t, end_t, SP) % Define output/inputs 
% function [totalspk] = countspk(start_t, end_t, SP) % Define output/inputs 
%
% Function to obtain number of spikes between start_t and end_t for each
% neuron in SP
%
% Inputs:
% start_t: start time in seconds (single number)
% end_t: end time in seconds (single number)
% SP: spike times, loaded from AllSpikes.mat (npxl pipeline output)
%
% Outputs:
% totalspk: 1 x nCells vector of spike counts
%
% 2023-07-11 Maya & Matt

    %SP.t_sec = (0:(Rows(SP.t_uS)-1))/1e-6; %Convert microseconds to seconds
    totalspk = zeros(size(SP)); % Make empty variable same length of dataset

    for i = 1:length(SP) %For every timepoint in the spike time dataset
        this_spk = SP(i).t_uS; 
        this_spk_s = this_spk * 10^-6;

     output =  this_spk_s > start_t & this_spk_s < end_t; %find the spikes between these timepoints
     totalspk(i) = sum(output); % Add up the spikes within these timepoints
    end
totalspk = totalspk'; %Invert so is column, not row   

end

