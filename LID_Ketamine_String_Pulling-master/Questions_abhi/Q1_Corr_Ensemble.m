GIX = sum(AllQ_all)' >200 & TBL.Depth_uM > 3000 & categorical(TBL.Group) == '6OHDA_LID';
% GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
% GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
sum(GIX), sum(GIXleft), sum(GIXright)


% Corr between the Baseline LDOPA and Ket periods

[c1, p1] = corr(ALL_R(:,1), ALL_R(:,2), 'rows', 'complete');

[c2, p2] = corr(ALL_R(:,1), ALL_R(:,3), 'rows', 'complete');

corrb = corr(AllQ_b(:,GIXright));
t = triu(ones(size(corrb)),1)==1;
corrb = corrb(t);
corrbl = corr(AllQ_b(:,GIXleft));
t = triu(ones(size(corrbl)),1)==1;
corrbl = corrbl(t);

corrp1 = corr(AllQ_p1(:,GIXright));
t = triu(ones(size(corrp1)),1)==1;
corrp1 = corrp1(t);
corrp1l = corr(AllQ_p1(:,GIXleft));
t = triu(ones(size(corrp1l)),1)==1;
corrp1l = corrp1l(t);

corrp2 = corr(AllQ_p2(:,GIXright));
t = triu(ones(size(corrp2)),1)==1;
corrp2 = corrp2(t);
corrp2l = corr(AllQ_p2(:,GIXleft));
t = triu(ones(size(corrp2l)),1)==1;
corrp2l = corrp2l(t);

[cp1, p1] = corr(corrb, corrp1, 'rows', 'complete');
[cp2, p2] = corr(corrb, corrp2, 'rows', 'complete');
[cp2p1, p2p1] = corr(corrp1, corrp2, 'rows', 'complete');

[cp1l, p1l] = corr(corrbl, corrp1l, 'rows', 'complete');
[cp2l, p2l] = corr(corrbl, corrp2l, 'rows', 'complete');


