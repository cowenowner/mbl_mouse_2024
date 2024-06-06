function OUT = PHO_acquire_data(ni_dev)
% Send out pulse train...
d = daq("ni");
d.Rate = 500;
ch1 = addinput(d,"Dev1","ai0","Voltage");
ch2 = addoutput(d,"Dev1","ao0","Voltage");
% output signal
% outputData = (linspace(-1,1,5000)');
% preload(d,outputData);
% start(d);
% Read data
DAQ_1 = read(d,seconds(1));
% Plot data
plot(DAQ_1.Time, DAQ_1.Variables)
xlabel("Time");ylabel('V')

