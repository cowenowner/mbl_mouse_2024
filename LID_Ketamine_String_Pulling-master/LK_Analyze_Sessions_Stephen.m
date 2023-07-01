%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads the master session database and allows
% you to select the subset of sets to analyze and the analysis
% function.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
global DIRS % This can be set up in the starup.m. See this file

ALLSES = LK_Load_SessionList_File(DIRS.SessionList_File);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter the analysis code from the Questions folder.
% AND select the appropriate sessions to analyze.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ana_type = 4;
ana_type = 3.1;
switch ana_type
    case 1
        ana_function = @Q1_What_does_ketamine_do_firing_rate; %Ketamine or ketamine = saline but not with Ldopa as that's different.
        GIX = ALLSES.Sorted & ALLSES.Session > 0 & ALLSES.Ketamine == 1 & ALLSES.LDOPA == 0; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.
    case 2
        ana_function = @Q5_What_does_LDOPA_do_firing_rate;
        GIX = ALLSES.Sorted & ALLSES.Session > 0 & ALLSES.LDOPA == 1; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.
    case 3
        ana_function = @Q1_What_does_ketamine_do_firing_rate;
        GIX = ALLSES.Sorted & ALLSES.Session > 0 & ALLSES.LDOPA == 1 & ALLSES.Ketamine == 0; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.
    case 3.1 % 80 Hz and LID
        ana_function = @Q3_Does_LDOPA_create_80Hz;
        GIX =  ALLSES.RatType == '6ODHA_LID' & ALLSES.LDOPA == 1; %
    case 4
        ana_function = @Q4_Do_spikes_phase_lock_to_LFP;
        GIX =  ALLSES.Sorted & ALLSES.LDOPA; %  
    case 5
        ana_function = @LK_Create_LFP_ratings_file;
        GIX =  ALLSES.Sorted & (ALLSES.Ketamine | ALLSES.LDOPA); %
    case 6
        ana_function = @LK_Determine_good_pull_bouts;
        GIX = ALLSES.StringPulling & ALLSES.Session > 0 ; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.
    case 7 % String
        ana_function = @Q1_Do_neurons_correlate_with_string_pulling_single_unit;
        GIX = ALLSES.Sorted & ALLSES.StringPulling & ALLSES.DeepLabCutFile & ALLSES.Session > 0 ; % Most sessions > 18 are current messed up due to errors with the naming conventions or missing files.
    otherwise
        error('wtf')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the chosen function on ALL of these sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(ALLSES(GIX,:))
LK_Session_Iterator(ana_function, ALLSES(GIX,:))
