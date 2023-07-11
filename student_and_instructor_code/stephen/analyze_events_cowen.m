% This should have the tcat.nidq.bin file.
data_dir = 'E:\NSB_Mouse\23242\HC101_23242_g0';
% This should be the CatGT command line or extracting the event times for
% an INVERTED pulse (like we have for the manually scored events).
% The -xia is inverted.

%For this data set where input 6 is the inverted hand-made inveted ttl pulse for events.
% Using SpikeGLX_NSIM I found that the INVERTED threshold would be around
% 1.3V. C:\Data\NSandB_course2023\HC101_23242_g0
%%
% Here is the command for the hand-coded pulses with the inverted by-hand pulses.
full_cmd = 'C:\CatGT-win\CatGT.exe -dir=C:\Data\NSandB_course2023 -run=HC101_23242 -g=0 -t=0 -ni -xia=0,0,6,-1.3,-1.0,0';

[status,cmdout] = system(full_cmd,'-echo');

% load and verify the times.
Event_times_sec = load('C:\Data\NSandB_course2023\HC101_23242_g0\HC101_23242_g0_tcat.nidq.xia_6_0.txt');

figure
plot(Event_times_sec,zeros(size(Event_times_sec)),'+')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now let's try a dataset with those very small inverted ttl pulses.
% D:\NSB_Mouse\23242\02_7_7_23\HC_2_Neuro2.0_g0\HC_2_Neuro2.0_g0_imec0 HC101_23242_g0_t0.nidq
% TRICK: Since these are negative going potentials, the threshold values
% should be NEGATIVE EVEN THOUGH the -xia flips the voltage.
full_cmd = 'C:\CatGT-win\CatGT.exe -dir=D:\NSB_Mouse\23242\02_7_7_23 -run=HC_2_Neuro2.0 -g=0 -t=0 -ni -inarow=52 -xia=0,0,6,-0.5,-0.3,0';


[status,cmdout] = system(full_cmd,'-echo');

% load and verify the times.
% C:\Data\NSandB_course2023\HC101_23242_g0
Event_times_sec = load('D:\NSB_Mouse\23242\02_7_7_23\HC_2_Neuro2.0_g0\HC_2_Neuro2.0_g0_tcat.nidq.xia_6_0.txt');

figure
plot(Event_times_sec,zeros(size(Event_times_sec)),'+')

%%

% Load the AI data...
if 0 % for testing
    obj = SGLX_Class;
    NIO.meta = obj.ReadMeta('HC_2_Neuro2.0_g0_t0.nidq.meta','D:\NSB_Mouse\23242\02_7_7_23\HC_2_Neuro2.0_g0');
    n_samples = 400e6; % obscene number to load all records.
    start_rec = 0;
    [NIO.data] = obj.ReadBinVolts(start_rec,n_samples,NIO.meta ,'HC_2_Neuro2.0_g0_t0.nidq.bin','D:\NSB_Mouse\23242\02_7_7_23\HC_2_Neuro2.0_g0'); %(samp0, nSamp, meta, binName, path)
    figure
    plot(NIO.data(5,:))
end