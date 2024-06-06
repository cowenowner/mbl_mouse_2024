function [SP, INFO] = PHY_read_spikes(data_dir,sFreq, groups_to_load)
% Reads spikes from the output of Phy. Stores them in a structure SP
%
% COWEN 2023 - Replaced loading the synced_spike_seconds.npy with
% spike_seconds.npy as the syncing was not required - mistake as this is
% already being done with the event files.
%
SP = []; % empty in case no spikes are found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    data_dir = pwd;
end
if nargin < 3
    groups_to_load = {'good'}; % you can also include mua.
end

if ~exist(fullfile(data_dir,'probe_depths.mat'),'file')
    prompt = {'Enter max depth of electrode array (tip) in millimeters(mm):'};
    dlgtitle = 'Array Depth in MM';
    dims = [1 35];
    definput = {''};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    if isempty(answer{1})
        error('you must enter a depth for the array in mm')
    end
    tip_depth_mm = str2double(answer{1});
    save(fullfile(data_dir,'probe_depths.mat'),'tip_depth_mm')
else
    load(fullfile(data_dir,'probe_depths.mat'),'tip_depth_mm');
end
if ~exist(fullfile(data_dir,'spike_seconds.npy'),'file')
    error('No spike_seconds.npy file found')
end

C = readNPY(fullfile(data_dir,'spike_clusters.npy'));
AMP = readNPY(fullfile(data_dir,'amplitudes.npy'));
% The 'synced_spike_seconds.npy' is produced by running
% NPXL_Sync_With_TPrime. You sadly need to do this to guarantee that the
% spike times are truly aligned.
T_usec = readNPY(fullfile(data_dir,'spike_seconds.npy'));
T_usec = T_usec*1e6;
if any(diff(T_usec)<0)
    disp('WARNING: spike_seconds is not in the correct order. Will attempt to correct.')
    sum(diff(T_usec)<0)
    % Correcting: This should fix all of the issues.
    [T_usec_sorted, six] = sort(T_usec); 
    C = C(six);
    AMP = AMP(six);
end

TPLT = readNPY(fullfile(data_dir,'templates.npy'));
TPLTind = readNPY(fullfile(data_dir,'templates_ind.npy')); % I don't think I need this, but need to be careful.
% INFO.Kilosort_REZ = load(fullfile(data_dir,'rez.mat')); % this is HUGE.
% Not sure why. It was never that big in earlier versions.
INFO.sFreq = sFreq;
INFO.data_dir = data_dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_usec = 1e6*double(T)/sFreq;
% Why are there fewer waveforms than spike times? Sould be the same!
% probably why 'subset' is in the name. How to get the recID for each one?
% It is  in the ch column in CI I believe.
% This file is not in the output - I did not run extract-waveforms so that
% is probably why.
WV = readNPY(fullfile(data_dir,'_phy_spikes_subset.waveforms.npy'));
% SC = readNPY(fullfile(data_dir,'_phy_spikes_subset.channels.npy'));
WV_recID = readNPY(fullfile(data_dir,'_phy_spikes_subset.spikes.npy'));
T_WV_usec = T_usec(WV_recID+1);

% CG = readtable(fullfile(data_dir,'cluster_group.tsv'),'FileType','text', 'Delimiter' ,'tab');
CI = readtable(fullfile(data_dir,'cluster_info.tsv'),'FileType','text', 'Delimiter' ,'tab');
CIgroup = categorical(CI.group);
cnt = 1;
for iG = 1:length(groups_to_load)
    gix = find(CIgroup == groups_to_load{iG}); % Rows in the cluster table from phy
    for ii = 1:length(gix)
        row_ix = gix(ii); % A particular row in the table.
        cid = CI.cluster_id(row_ix);
        SPK_IX = C == cid;

        t_uS = T_usec(SPK_IX);
        amp = single(AMP(SPK_IX));
        
        if any(diff(t_uS) <=0)
            disp('WARNING: out of order or duplicate timestamps')
            % find these timestamps in SPK_IX and remove 
            [t_uS,ia] = unique(t_uS);
            amp = amp(ia);
        end

        SP(cnt).t_uS = t_uS;
        SP(cnt).cluster_id = cid;
        SP(cnt).Amplitude = CI.Amplitude(row_ix);
        SP(cnt).KSLabel = CI.KSLabel{row_ix};
        SP(cnt).PHYLabel = CI.group{row_ix};
        SP(cnt).amp = CI.amp(row_ix);
        SP(cnt).amp_all = amp;
        SP(cnt).template_index0 = CI.ch(row_ix); % this is the index to the template for this cluster (I am 87% certain of this).
        SP(cnt).depth_on_probe = CI.depth(row_ix);
        SP(cnt).depth_of_probe_tip_uM = tip_depth_mm*1000;
        SP(cnt).neuropixels_depth_uM = tip_depth_mm*1000 - CI.depth(row_ix); % THis is only appropriate for neuropixels. 1000 converts depth from uM to mm.
        SP(cnt).fr = CI.fr(row_ix);
        SP(cnt).n_spikes = CI.n_spikes(row_ix);
        SP(cnt).fname = fullfile(data_dir,'spike_clusters.npy'); % legacy
        % Template
        ix = find(TPLTind(1,:) == SP(cnt).template_index0);
        SP(cnt).template = single(squeeze(TPLT(:,:,ix))); % not sure
        % WV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Since the WV are subsampled - find the overlapping waveforms with the
        % reciDs.
        %         WV_ix = intersect(find(SPK_IX)-1,find(WV_recID)); % need -1 as recID starts at 0
        [~,WV_ix] = intersect(T_WV_usec,t_uS); % need -1 as recID starts at 0

        SP(cnt).WV.mWV = single(squeeze(mean(WV(WV_ix,:,:),1,'omitnan')));
        SP(cnt).WV.sWV = single(squeeze(std(WV(WV_ix,:,:),[],1,'omitnan')));
        SP(cnt).dir = data_dir;
        SP(cnt).Cell_type = [];
        SP(cnt).Quality = [];
        SP(cnt).Notes = [];
        SP(cnt).mfile = mfilename;
        SP(cnt).CQ = []; % cluster quality (TODO). Cluster separation from mclust was what it used to be.

        INFO.t_uS{cnt} =t_uS; % store the same thing as a cell array which can be convenient for plotting.
        cnt = cnt + 1;
    end
end

if nargout ==0
    T = Restrict(INFO.t_uS, 0, 20e6);
    for ii = 1:length(T)
        T{ii} = T{ii}/1e6;
    end
    figure
    plot_raster(T)
    axis tight
    ylabel('Neuron ID')
    xlabel('sec')
    pubify_figure_axis;

    %% Corr
    Q = Bin_ts_array(T,0.05);
    R = corr(Q);

    figure
    Plot_R_matrix(R);
end