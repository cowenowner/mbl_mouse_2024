function OUT = WV_long_meta_Abhi()
OUT = [];
% defining bin sizes and things
bin_size_ac_ms = 2;
n_ac_bins = 200;
neuron_quality_threshold = 2;

[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;
OUT.SP = SP;
% catch
%     disp('Failed to load files on first try for some reason. Box?')
%     pause(10)
%     [GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;
% end
OUT.n_ac_bins = n_ac_bins;
OUT.bin_size_ac_ms = bin_size_ac_ms;

OUT.DEPTHS = DEPTHS; % Send this info out of the function for meta analysis.
OUT.META = META;
OUT.EVT = EVT;

load('AllSpikes_longWV_butter.mat')
cnt = 1;
for iS = 1:length(SP)
    if SP(iS).Quality > neuron_quality_threshold
        SPtmp(cnt) = SP(iS);
        cnt = cnt + 1;
    end
end
SP = SPtmp;
OUT.SP_WV_long = SP;

%% Extract the relevant times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_start_end_uS = [];

if any(E.EventID == 'KetInjectionStart') && any(E.EventID == 'SalineInjectionStart')
        e_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
        e_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-55 -35; -25 -5; 2 30; 60 80];
        big_peri_event_min = [-55 80];

elseif any(E.EventID == 'LDOPAInjectionStart')
        e_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
        e_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 25 55; 65 90; 120 140];
        big_peri_event_min = [-25 140];
end

OUT.intervals_around_evt_min = intervals_around_evt_min;
TIMES.EventStartEndUsec(1,1) = e_start_end_uS(1);
TIMES.EventStartEndUsec(1,2) = e_start_end_uS(2);
TIMES.BaselineUsec(1,1) = TIMES.EventStartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.BaselineUsec(1,2) = TIMES.EventStartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.PostInjectionUsec(1,1) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(2,1)*60*1e6;
TIMES.PostInjectionUsec(1,2) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(2,2)*60*1e6;
TIMES.PostInjectionUsec(2,1) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(3,1)*60*1e6;
TIMES.PostInjectionUsec(2,2) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(3,2)*60*1e6;
TIMES.PostInjectionUsec(3,1) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(4,1)*60*1e6;
TIMES.PostInjectionUsec(3,2) = TIMES.EventStartEndUsec(2) + intervals_around_evt_min(4,2)*60*1e6;

TIMES.PeriKetamineUsec = [TIMES.EventStartEndUsec(1) + big_peri_event_min(1)*60*1e6 TIMES.EventStartEndUsec(1) + big_peri_event_min(2)*60*1e6 ];

OUT.TM = TIMES;

%Get autocorr
[AC,x] = AutoCorrArray_Abhi(TS,bin_size_ac_ms*1000,n_ac_bins*1000,[TIMES.BaselineUsec;TIMES.PostInjectionUsec]);
% ACd = AC -  AC(:,:,1);
% cm = lines(5);
OUT.AC_base_post1_post2_post3 = AC;
OUT.AC_x_ms = x/1000;
