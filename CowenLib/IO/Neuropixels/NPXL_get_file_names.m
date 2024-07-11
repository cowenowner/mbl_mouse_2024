function [D] = NPXL_get_file_names(data_folder,npxl_top_dir_name,kilosort_dir_name )
% function [D] = NPXL_get_file_names(data_folder,npxl_top_dir_name,kilosort_dir_name )
%
% Figure out the names and folders of important files in the Neuropixels
% data directories. You should probably tcat things as this assumes tcat
% files or at least gives priority to these files.
%
% Cowen 2024
D.ni_data_dir = fullfile(data_folder,npxl_top_dir_name);
D.ap_data_dir = fullfile(data_folder,npxl_top_dir_name,[npxl_top_dir_name '_imec0']);
D.kilosort_dir = fullfile(D.ap_data_dir,kilosort_dir_name);

d1 = dir(fullfile(D.ni_data_dir,'*tcat.nidq.bin')); %PhotoPixelsStrobe_g0_t0.nidq.bin
d2 = dir(fullfile(D.ni_data_dir,'*.nidq.bin')); %PhotoPixelsStrobe_g0_t0.nidq.bin
if ~isempty(d1)
    D.nidq_bin_file_path = fullfile(D.ni_data_dir,d1(1).name);
else
    D.nidq_bin_file_path = fullfile(D.ni_data_dir,d2(1).name);
end

d = dir(fullfile(D.ap_data_dir,'*tcat*.lf.bin'));
if isempty(d)
    D.lfp_bin_file_path = [];
    D.lfp_fname = [];
    D.meta_fname_lf = [];
    disp('No LFP data - possible the .lf.bin file was not created from the .ap.bin file.')
else
    D.lfp_bin_file_path = fullfile(D.ap_data_dir,d(1).name);
    [path_name, tmp] = fileparts(D.lfp_bin_file_path);
    D.lfp_fname = [tmp '.bin'];
    D.meta_fname_lf = strrep(D.lfp_fname,'.lf.bin','.lf.meta');
end

d = dir(fullfile(D.ni_data_dir,'synced*tcat.nidq*.txt'));
if isempty(d)
    d = dir(fullfile(D.ni_data_dir,'*tcat.nidq*.txt'));
end
if isempty(d)
    d = dir(fullfile(D.ni_data_dir,'*.nidq*.txt'));
    % 
    disp('No tcat files found. Using the original data.')
end

cnt = 1;
D.event_files = [];
for ii = 1:length(d)
    if d(ii).bytes > 0
        D.event_files{cnt} = fullfile(d(ii).folder,d(ii).name);
        cnt = cnt + 1;
    end
end
