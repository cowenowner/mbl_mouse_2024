function SP = PHY_Process_Intan_Sorted_Neurons_to_AllSpikes(phy_dir,GET_USER_INPUT)
%
% After you are done spike sorting in Phy, run this function and it will create a nice AllSpikes.mat files that have all of the
% relevant meta data like waveform, etc... Nice matlab files so no need for
% silly readers for specialized binary files. Also prompts for cluster
% quality.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2021, 2022 - adapted to neuropixels.
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
SP = PHY_read_spikes_intan(phy_dir,sFreq);

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
        saveas(gcf,fullfile(phy_dir,[nm '_clust_summary.png']))
        prompt = {'NUMERIC Quality Score:1=worst,5=best','Est. Neuron Type (PYR,INT)','Amp cutoff:', 'Notes:'};
        dlgtitle = ['Rate Cluster ' nm];
        dims = [1 35];
        definput = {'','','',''};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        SP(iF).Quality = [];
        SP(iF).Cell_Type = [];
        SP(iF).Amp_cutoff = [];
        %     TODO MAYBE:    SP(iF).ISI_cutoff_ms = 1.2; % get rid if spikes that are super close to one another as this is likely due to the spike sorting error.

        SP(iF).removed_using_amp_threshold = 0;
        SP(iF).Notes = [];

        if ~isempty(answer)
            SP(iF).Quality = str2double(answer{1});
            SP(iF).Cell_Type = answer{2};
            SP(iF).Amp_cutoff = str2double(answer{3});
            SP(iF).Notes = answer{4};
        end
        % apply an amplitude cutuff if desired to get rid of artifact super
        % high amplitude spikes if they exist.
        if ~(isempty(SP(iF).Amp_cutoff) || SP(iF).Amp_cutoff == 0 || isnan(SP(iF).Amp_cutoff))
            figure(20)
            clf
            plot(SP(iF).t_uS/1e6,SP(iF).amp_all,'.')
            hold on

            GIX = SP(iF).amp_all < SP(iF).Amp_cutoff;
            SP(iF).amp_all =  SP(iF).amp_all(GIX);
            SP(iF).t_uS =  SP(iF).t_uS(GIX);
            SP(iF).n_spikes =  sum(GIX);
            SP(iF).amp =  mean(SP(iF).amp_all);
            SP(iF).removed_using_amp_threshold = sum(GIX);
            
            plot(SP(iF).t_uS/1e6,SP(iF).amp_all,'o')
            xlabel('sec')

            fprintf('n spik < 1.2 ms ISI %d \n',sum(diff(SP(iF).t_uS/1000)<1.2))
        end
    end
end
save(fullfile(phy_dir,'AllSpikes.mat'),'SP')
disp('Saved as AllSpikes.mat')
if nargout == 0
    PHY_Check_Clusters_From_AllSpikes(SP)
end
