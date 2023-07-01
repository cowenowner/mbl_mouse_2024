DatFile1 = INTAN_Read_DAT_file('amp-A-042.dat');
Reduce_d1 = DatFile1(1:600000);
F = designfilt('bandpassiir','FilterOrder',12, ...
                    'HalfPowerFrequency1',250,'HalfPowerFrequency2',6000, ...
                    'SampleRate',30000);
 Filt_lfp = filtfilt(F,Reduce_d1);
 figure;plot(Filt_lfp(1:300000))