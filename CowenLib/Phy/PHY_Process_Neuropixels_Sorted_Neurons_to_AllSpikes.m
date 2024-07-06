function SP = PHY_Process_Neuropixels_Sorted_Neurons_to_AllSpikes(phy_dir,GET_USER_INPUT)
% function SP = PHY_Process_Neuropixels_Sorted_Neurons_to_AllSpikes(phy_dir,GET_USER_INPUT)
%
% After you are done spike sorting in Phy, run this function and it will create a nice AllSpikes.mat files that have all of the
% relevant meta data like waveform, etc... Also prompts for cluster
% quality.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    phy_dir = pwd;
end
if nargin < 2
    GET_USER_INPUT = true; % don't do this if you just want a quick AllSpikes file for testing.
end
% REMOVE_SHORT_ISIs = true;
ISI_thresh_ms = 1;
meta_file = dir(fullfile(phy_dir,'*tcat.*.ap.meta'));
if isempty(meta_file)
    % assume the meta is in the directory above.
    meta_file = dir(fullfile(phy_dir,'..','*tcat.*.ap.meta'));
end

if length(meta_file) >1 
    phy_dir
    error('MORE THAN ONE META FILE!')
end
% Assume an rhd file and load the spikes from the Phy output.
% Load meta data.
obj = SGLX_Class;
meta = obj.ReadMeta(meta_file(1).name,meta_file(1).folder);
sFreq = str2double(meta.imSampRate); % MAKE SURE THAT THIS IS THE TRUE RATE!!!!
if sFreq == 30000
    disp('WARNING: The sample rate in your ap.meta file is exactly 30000. THIS IS NOT THE TRUE RATE. Run the calibration tool in SPikeGLX and update the .meta file yourself. https://billkarsh.github.io/SpikeGLX/help/syncEdges/Sync_edges/')
end

SP = PHY_read_spikes(phy_dir,sFreq,{'good'}); % This is from cluster_info.txv. The last col is the group which is the user-specified rating (green rows in phy).
disp(['Found ' num2str(length(SP)) ' possible cells'])

for iF = 1:length(SP)
    GIX = diff(SP(iF).t_uS/1000) > ISI_thresh_ms;
    fprintf('%d %2.2f percent of spikes have short ISIs\n',iF, 100*(1-(sum(GIX))/length(SP(iF).t_uS)));
    SP(iF).IX_of_good_spikes = GIX;
    
    nm = sprintf('cid_%d',SP(iF).cluster_id);
    tmp = SP(iF);
    PHY_Check_Clusters_From_AllSpikes(tmp)
    saveas(gcf,fullfile(phy_dir,[nm '_clust_summary.png']))

    if GET_USER_INPUT
        prompt = {'NUMERIC Quality Score:1=worst,5=best','Est. Neuron Type (PYR,INT)','Notes:'};
        dlgtitle = [num2str(iF) ') Rate Cluster ' nm];
        dims = [1 35];
        definput = {'2','',''};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        SP(iF).Quality = str2double(answer{1});
        SP(iF).Cell_Type = answer{2};
        SP(iF).Notes = answer{3};
    end
    fprintf('%d ',iF)
end
save(fullfile(phy_dir,'AllSpikes.mat'),'SP')
disp('Saved as AllSpikes.mat')
if nargout == 0
    PHY_Check_Clusters_From_AllSpikes(SP)
end
