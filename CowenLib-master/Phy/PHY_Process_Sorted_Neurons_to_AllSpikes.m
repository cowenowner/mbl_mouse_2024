function SP = PHY_Process_Sorted_Neurons_to_AllSpikes(phy_dir,GET_USER_INPUT)
%
% After you are done spike sorting in Phy, run this function and it will create a nice AllSpikes.mat files that have all of the
% relevant meta data like waveform, etc... Nice matlab files so no need for
% silly readers for specialized binary files. Also prompts for cluster
% quality.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    phy_dir = pwd;
end
if nargin < 2
    GET_USER_INPUT = true;
end

% Assume an rhd file and load the spikes from the Phy output.
RHD = read_Intan_RHD2000_file_cowen(phy_dir,'info.rhd');
sFreq = RHD.frequency_parameters.amplifier_sample_rate;
SP = PHY_read_spikes(phy_dir,sFreq);

for iF = 1:length(SP)
    % Spikes can't be within 1ms of each other. Would indicate overlapping
    % spikes and I would not trust the spike sorting for such spikes given waveform superposition. Get
    % rid of them.
    %     GIX = diff(SPK(iF).t_uS/1000) > 1;
    %     SPK(iF).t_uS = SPK(iF).t_uS(GIX);
    nm = sprintf('cid_%d',SP(iF).cluster_id);
    if GET_USER_INPUT
        tmp = SP(iF);
        PHY_Check_Clusters_From_AllSpikes(tmp)
        prompt = {'NUMERIC Quality Score:1=worst,5=best','Est. Neuron Type (PYR,INT)','Notes:'};
        dlgtitle = ['Rate Cluster ' nm];
        dims = [1 35];
        definput = {'','',''};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        SP(iF).Quality = str2double(answer{1});
        SP(iF).Cell_Type = answer{2};
        SP(iF).Notes = answer{3};
    end
end
save(fullfile(phy_dir,'AllSpikes.mat'),'SP')
disp('Saved as AllSpikes.mat')
if nargout == 0
    PHY_Check_Clusters_From_AllSpikes(SP)
end
