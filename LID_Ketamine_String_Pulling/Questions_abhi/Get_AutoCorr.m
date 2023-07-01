function OUT = Get_AutoCorr()

OUT = [];
bin_size_ac_ms = 2;
acorr_width_ms = 200; 

[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things;

OUT.SP = SP;
OUT.Event_times = E;

[AC,x] = AutoCorrArray(TS,bin_size_ac_ms*1000,acorr_width_ms*1000);

OUT.AC = AC;
OUT.AC_x_ms = x/1000;