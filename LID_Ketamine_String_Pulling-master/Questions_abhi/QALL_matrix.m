function OUT = QALL_matrix()

OUT = [];
binsize_Q_ms = 20;
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;

inj1_start_end_uS = [];
inj2_start_end_uS = [];

if any(E.EventID == 'KetInjectionStart') && any(E.EventID == 'SalineInjectionStart')
        inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
        inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
        inj1_start_end_uS(1) = E.MinFromStart(E.EventID == 'SalineInjectionStart')*60*1e6;
        inj1_start_end_uS(2) = E.MinFromStart(E.EventID == 'SalineInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 2 22; 60 80];
        big_peri_event_min = [-25 80];
        OUT.inj1 = 'Sal';
        OUT.inj2 = 'Ket';

elseif any(E.EventID == 'LDOPAInjectionStart')
        inj1_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
        inj1_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        if any(E.EventID == 'SalineInjectionStart')
            inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'SalineInjectionStart')*60*1e6;
            inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'SalineInjectionEnd')*60*1e6;
            OUT.inj2 = 'Sal';
        elseif any(E.EventID == 'KetInjectionStart')
            inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
            inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
            OUT.inj2 = 'Ket';
        end
        intervals_around_evt_min =  [-25 -5; 2 22; 60 80; 5 25];
        big_peri_event_min = [-25 180];
        OUT.inj1 = 'Ldo';
end

% return
OUT.intervals_around_evt_min = intervals_around_evt_min;
TIMES.Event1StartEndUsec(1,1) = inj1_start_end_uS(1);
TIMES.Event1StartEndUsec(1,2) = inj1_start_end_uS(2);
TIMES.Event2StartEndUsec(1,1) = inj2_start_end_uS(1);
TIMES.Event2StartEndUsec(1,2) = inj2_start_end_uS(2);
TIMES.BaselineUsec(1,1) = TIMES.Event1StartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.BaselineUsec(1,2) = TIMES.Event1StartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.PostInjectionUsec(1,1) = TIMES.Event2StartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.PostInjectionUsec(1,2) = TIMES.Event2StartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.PostInjectionUsec(2,1) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(2,1)*60*1e6;
TIMES.PostInjectionUsec(2,2) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(2,2)*60*1e6;
TIMES.PostInjectionUsec(3,1) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(3,1)*60*1e6;
TIMES.PostInjectionUsec(3,2) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(3,2)*60*1e6;


TIMES.PeriKetamineUsec = [TIMES.Event1StartEndUsec(1) + big_peri_event_min(1)*60*1e6 TIMES.Event1StartEndUsec(1) + big_peri_event_min(2)*60*1e6 ];

edges = TIMES.PeriKetamineUsec(1):binsize_Q_ms*1000:TIMES.PeriKetamineUsec(2);
[Qall,Qall_x_uS] = Bin_ts_array(TS, edges);
Qall = int16(Qall);

OUT.Qall = Qall;
OUT.Qall_x_uS = Qall_x_uS(:,1) + binsize_Q_ms*1000/2;
