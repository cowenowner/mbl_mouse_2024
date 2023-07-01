%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads the master session database and allows 
% you to select the subset of sets to analyze and the analysis 
% function.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
global DIRS
% Set up the locations for all of the important directories and files.
% Do not create new DIRS varibles unless we all agree upon them
% DIRS.Data_Dir = 'C:\Users\Stephen Cowen\Box\Cowen Laboratory\Data\LID_Ketamine_Single_Unit_R56';

DIRS.Data_Dir = 'E:\META_analysis';
DIRS.Box_Dir = 'C:\Users\Stephen Cowen\Box\Cowen Laboratory';
DIRS.Analysis_Dir = 'E:\Temp_6.6.23';
DIRS.SessionList_File = 'C:\Users\Stephen Cowen\Box\Cowen Laboratory\Data\LID_Ketamine_Single_Unit_R56\SessionInfo.xlsx';
DIRS.LFP_Dir = '';
DIRS.Video_Dir = '';


ALLSES = LK_Load_SessionList_File(DIRS.SessionList_File);
%GP = LK_Globals;
% You can override the assumed folders created by LK_Globals.
% GP.Data_Dir = 'C:\Users\abhi2\Box Sync\Cowen Laboratory\Data\LID_Ketamine_Single_Unit_R56'; % abhi laptop
% %GP.Data_Dir = 'C:\Users\Stephen Cowen\Box Sync\Cowen Laboratory\Data\LID_Ketamine_Single_Unit_R56'; % The top directory where all of the sessions live..
% GP.Analysis_Dir = 'C:\Temp\TempAnaResults'; % A full path to a folder (must exist) where all of the results go.
% GP.SessionList_File = 'C:\Users\abhi2\Box Sync\Cowen Laboratory\!Projects\LID_Ketamine_Single_Unit_R56\Meta_Data\SessionInfo.xlsx'; % abhi laptop
% %GP.SessionList_File = 'C:\Users\Stephen Cowen\Box Sync\Cowen Laboratory\!Projects\LID_Ketamine_Single_Unit_R56\Meta_Data\SessionInfo.xlsx';

%ALLSES = LK_Load_SessionList_File(GP.SessionList_File);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter the analysis code from the Questions folder.
% AND select the appropriate sessions to analyze.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ana_function = @ Q19_Is_there_cross_freq_coupling_LID_ket_Abhi;
GIX = ALLSES.Sorted; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ana_function = @Q1_sal_vs_ketamine_single_unit;
% GIX = ALLSES.Sorted & ALLSES.Saline & ALLSES.Ketamine;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ana_function = @Q1_What_does_ketamine_do_firing_rate;
% GIX = ALLSES.Sorted & ALLSES.Ketamine & ALLSES.Session > 30 & ALLSES.Rat == 315 ; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3%%%%%%%%%%%%%%
% ana_function = @Q2_Do_neurons_correlate_with_string_pulling;
% %  GIX = ALLSES.StringPulling;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the chosen function on ALL of these sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running on...')
ALLSES(GIX,:)
LK_Session_Iterator(ana_function, ALLSES(GIX,:))
