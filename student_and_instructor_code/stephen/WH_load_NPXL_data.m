function [D] = WH_load_NPXL_data(data_folder,npxl_top_dir_name,kilosort_dir_name )
% Figure out the names and folders of important files.
% Cowen 2024
D.ni_data_dir = fullfile(data_folder,npxl_top_dir_name);
D.ap_data_dir = fullfile(data_folder,npxl_top_dir_name,[npxl_top_dir_name '_imec0']);
D.kilosort_dir = fullfile(D.ap_data_dir,kilosort_dir_name);

d = dir(fullfile(D.ap_data_dir,'*tcat*.lf.bin'));
D.lfp_bin_file_path = fullfile(D.ap_data_dir,d(1).name);
[path_name, tmp] = fileparts(D.lfp_bin_file_path);
D.lfp_fname = [tmp '.bin'];
D.meta_fname_lf = strrep(D.lfp_fname,'.lf.bin','.lf.meta');
d = dir(fullfile(D.ni_data_dir,'synced*tcat.nidq*.txt'));
cnt = 1;
D.event_files = [];
for ii = 1:length(d)
    if d(ii).bytes > 0
        D.event_files{cnt} = fullfile(d(ii).folder,d(ii).name);
        cnt = cnt + 1;
    end
end
