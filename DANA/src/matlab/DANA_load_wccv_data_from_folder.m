function [BLOCK_INFO, CV_DATA] = DANA_load_wccv_data_from_folder(WCCV_data_folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the CV data and stimulation times from the excel sheet.
%
% Data expected to be in a folder where each text file is a 2 col tab
% separated matrix of scan times in seconds (col 1) and CV values (col 2).
% The filename is <Fq>Hz_LV_<LV val>_<dist name>_<trial_num>_<stuff>.txt
%
% It also expects a xlsx file that is 3 col that has the list of each
% file listed above (col 1), block number (col 2), MM:SS (col 3)
%
% 20Hz_LV_1.01_gamma_	1	12:20
% 10Hz_LV_0.39_gamma_	1	12:25
% 10Hz_LV_0.00_WCCV_	1	12:30
% 20Hz_LV_0.00_WCCV_	1	12:40
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cowen 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = dir(fullfile(WCCV_data_folder,'*.xlsx'));
T = readtable(fullfile(WCCV_data_folder,d(1).name));
% Remove spaces from file names if there are any...
dtxt = dir(fullfile(WCCV_data_folder,'*.txt'));
for iD = 1:length(dtxt)
    if length(dtxt(iD).name) ~= length(strrep(dtxt(iD).name,' ',''))
        newname = strrep(dtxt(iD).name,' ','');
        disp(['Moving ' newname ])
        movefile(fullfile(dtxt(iD).folder, dtxt(iD).name), fullfile(dtxt(iD).folder, newname))
    end
end

CV_DATA = [];
BLOCK_INFO = [];

last_val = 0;

for iR = 1:Rows(T)
    % the filename prefix from the xlsx file...
    name = char(T{iR,1});
    % Parse the Hz and Lv from the prefix
    ix = strfind(name,'_');
    HZ = str2double(name(1:(ix(1)-3)));
    LV = str2double(name((ix(2)+1):(ix(3)-1)));
    % Get the associated CV file
    tmp = fullfile(WCCV_data_folder,[name num2str(T{iR,2}) '_*.txt']);
    tmp = strrep(tmp,' ','');
    dd = dir(tmp);
    tmp = readtable(fullfile(WCCV_data_folder,dd(1).name));
    % collect the data in a 2 col matlab matrix.
    new_data = [tmp.Var1 + last_val + tmp.Var1(2) tmp.Var2 ];
    % Append important related data as additional columns ...
    new_data = [new_data repmat([HZ LV],size(new_data,1),1)];

    CV_DATA = [CV_DATA; new_data ];

    last_val = new_data(end,1); % update the last timestamp so that we know when to start the times from the next block.

    % Store some meta-data for each block
    BLOCK_INFO(iR).fname_prefix = char(T{iR,1});
    BLOCK_INFO(iR).Hz = HZ;
    BLOCK_INFO(iR).LV = LV;
    BLOCK_INFO(iR).within_block_order = T{iR,2};
    BLOCK_INFO(iR).within_block_CV_data= [tmp.Var1 tmp.Var2];
    %     BLOCK_INFO(iR).block_start_end_sec = new_data([1 end],1)';
    BLOCK_INFO(iR).block_dur_sec = diff(new_data([1 end],1));
end
%
if nargout == 0
    figure
    plot(CV_DATA(:,1), CV_DATA(:,2))
    xlabel('sec')
    ylabel('[DA] units?')
end
