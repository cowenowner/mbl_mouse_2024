%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add appropriate paths including this directory and the CowenLib
% Directory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0, 'defaultTextInterpreter', 'none');
set(0, 'defaultLegendInterpreter', 'none');
set(0, 'defaultAxesTickLabelInterpreter', 'none');
set(0, 'defaultAxesFontName', 'Arial');
format short g

addpath(genpath(pwd))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cowenlibdir = fullfile(pwd,'..','CowenLib');

if exist(cowenlibdir,'dir')
    addpath(genpath(cowenlibdir))
else
    error('Could not find the CowenLib directory on your computer.')
end

dbstop if error
set(0, 'defaultFigureColormap', jet);

if contains(pwd,'Stephen')
    global DIRS
    % Set up the locations for all of the important directories and files.
    % Do not create new DIRS varibles unless we all agree upon them
    DIRS.Data_Dir = 'C:\Users\Stephen Cowen\Box\Cowen Laboratory\Data\LID_Ketamine_Single_Unit_R56';
    DIRS.Box_Dir = 'C:\Users\Stephen Cowen\Box\Cowen Laboratory';
    DIRS.Analysis_Dir = 'C:\Users\Stephen Cowen\Dropbox\Foldershare\Analysis_Results_Dropbox';
    DIRS.SessionList_File = 'C:\Users\Stephen Cowen\Box\Cowen Laboratory\!Projects\LID_Ketamine_Single_Unit_R56\Meta_Data\SessionInfo.xlsx';
    DIRS.LFP_Dir = 'Z:\Data\LID_Ketamine_Single_Unit_R56\LFP_Data';
    DIRS.Video_Dir = 'D:\R56_Single_Unit\Videos';
    disp (' Added DIRS global variable for Cowen')
end

disp('>>>> Added LID and CowenLib paths. <<<<')