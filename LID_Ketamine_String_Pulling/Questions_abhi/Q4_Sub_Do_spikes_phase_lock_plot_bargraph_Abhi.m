% bar plots of sig modulated neurons to phase and to power during peak 80
% Hz oscillation
IX_80_LDO = categorical(TBL.Condition) == 'LDO' & TBL.Interval == 3;
IX_50_KET = categorical(TBL.Condition) == 'KET' & TBL.Interval == 2;
x = categorical({'Peak 80 Hz','Peak 50 Hz'});
x = reordercats(x,{'Peak 80 Hz','Peak 50 Hz'});
y = [8 1; 16 5];
figure
bar(x,y,'stacked')

x1 = categorical({'Significant phase modulation','Significant power modulation'});
x1 = reordercats(x1,{'Significant phase modulation','Significant power modulation'});
y1 = [1 5; 8 16];
figure
bar(x1,y1)