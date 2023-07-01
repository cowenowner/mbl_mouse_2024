%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads the master session database and allows 
% you to select the subset of sets to analyze and the analysis 
% function.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GP = LK_Globals;
GP.Data_Dir = 'C:\Users\gowri\Box\Cowen Laboratory\Data\LID_Ketame_Single_Unit_R56'; 
GP.Analysis_Dir = 'C:\Users\gowri\OneDrive\Documents\CowenLab\Analysis'; % A full path to a folder (must exist) where all of the results go.
GP.SessionList_File = 'C:\Users\gowri\Box\Cowen Laboratory\!Projects\LID_Ketamine_Single_Unit_R56\Meta_Data\SessionInfo.xlsx'; 

ALLSES = LK_Load_SessionList_File(GP.SessionList_File);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter the analysis code from the Questions folder.
% AND select the appropriate sessions to analyze.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ana_function = @Q1_What_does_ketamine_do_firing_rate;
% GIX = ALLSES.Sorted & ALLSES.Ketamine & ALLSES.Session > 0 ; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ana_function = @Q2_Does_paw_track_corr_with_pos_imu;
GIX = ALLSES.Sorted & ALLSES.StringPulling & ALLSES.DeepLabCutFile & ALLSES.Session > 0 ; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.

% ana_function = @LK_Determine_good_pull_bouts;
% GIX = ALLSES.StringPulling & ALLSES.Session > 0 ; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the chosen function on ALL of these sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running on...')
ALLSES(GIX,:)
LK_Session_Iterator(ana_function, ALLSES(GIX,:),GP.Data_Dir,GP.Analysis_Dir)
