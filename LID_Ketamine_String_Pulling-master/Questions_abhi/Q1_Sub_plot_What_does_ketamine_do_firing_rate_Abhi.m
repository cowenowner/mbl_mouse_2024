%%% Modulation of neurons to ketmaine %%%
% Control M1 
IX_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control';
IX_Lesion_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Hemisphere == 'R' ;
IX_UnLesion_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Hemisphere == 'L';
% & categorical(TBL.Drugs) == 'LDOPA&Ketamine' 
IX_Lesion_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Hemisphere == 'R' ...
    & categorical(TBL.Drugs) == 'LDOPA&Ketamine';
IX_UnLesion_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Hemisphere == 'L'...
    & categorical(TBL.Drugs) == 'LDOPA&Ketamine';

IX_CS = categorical(TBL.BrainRegion) == 'Striatum' & categorical(TBL.Group) == 'Control';

% creating needed logicals
IX_LK = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine';
IX_SK = categorical(TBL_M1.Drugs) == 'Saline&Ketamine';
IX_LKLS = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' | categorical(TBL_M1.Drugs) == 'LDOPA&Saline';
IX_UnLes = categorical(TBL_M1.Hemisphere) == 'L';
IX_Les = categorical(TBL_M1.Hemisphere) == 'R';
IX_SK_Les = IX_SK & IX_Les;
IX_SK_UnLes = IX_SK & IX_UnLes;
IX_LK_Les = IX_SK & IX_Les;
IX_LK_UnLes = IX_SK & IX_UnLes;

% Looking at auto corr



% looking at LV before and after ketamine
figure
gscatter(TBL_M1.LocVar_base(IX_SK),TBL_M1.LocVar_post2(IX_SK,:),TBL_M1.NeuronType(IX_SK),clrs(1:3,:),'o+x',5)
figure
gscatter(TBL_M1.LocVar_base(IX_LK),TBL_M1.LocVar_post2(IX_LK,:),TBL_M1.NeuronType(IX_LK),clrs(1:3,:),'o+x',5)
xlabel('LV baseline')
ylabel('LV ketamine')
title('All LID rats all neurons LDOPA and Ketmaine session IN(1) PY(2) Regular spiking PY(3)')
pubify_figure_axis

figure
gscatter(TBL_M1.LocVar_base(IX_SK_UnLes),TBL_M1.LocVar_post2(IX_SK_UnLes),TBL_M1.NeuronType(IX_SK_UnLes),clrs(1:3,:),'o+x',5)
figure
gscatter(TBL_M1.LocVar_base(IX_LK_UnLes),TBL_M1.LocVar_post2(IX_LK_UnLes),TBL_M1.NeuronType(IX_LK_UnLes),clrs(1:3,:),'o+x',5)
xlabel('LV baseline')
ylabel('LV ketamine')
title('All LID rats UnLesioned hemisphere neurons Saline and Ketmaine session IN(1) PY(2) Regular spiking PY(3)')
pubify_figure_axis

% Looking at bursty cells LV > 1 during baseline 
IX_burst = IX_SK & TBL_M1.LocVar_base>1;
figure
gscatter(TBL_M1.LocVar_base(IX_burst),TBL_M1.LocVar_post2(IX_burst),TBL_M1.NeuronType(IX_burst),clrs(1:3,:),'o+x',5)

% Local Variance and depth
figure
scatter(TBL_M1.Depth_uM(IX_UnLes), TBL_M1.LocVar_base(IX_UnLes))
xlabel('Depth uM')
ylabel('LV')
title('All neurons from Lesioned hemisphere of LID rats LocVar during baseline for the depths in M1')
pubify_figure_axis

figure
scatter(TBL_M1.Depth_uM(IX_LKLS), TBL_M1.LocVar_base(IX_LKLS))
hold on
scatter(TBL_M1.Depth_uM(IX_LKLS), TBL_M1.LocVar_post1(IX_LKLS))
scatter(TBL_M1.Depth_uM(IX_LKLS), TBL_M1.LocVar_post2(IX_LKLS))

figure
scatter(TBL_M1.Depth_uM(IX_LK), TBL_M1.LocVar_base(IX_LK))
hold on
scatter(TBL_M1.Depth_uM(IX_LK), TBL_M1.LocVar_post1(IX_LK))

scatter(TBL_M1.Depth_uM(IX_LK), TBL_M1.LocVar_post2(IX_LK))
% hold on
% P = polyfit(TBL_M1.Depth_uM(IX_LK), TBL_M1.LocVar_post2(IX_LK), 1); % getting the linear fit
% Bfit = polyval(P, TBL_M1.Depth_uM(IX_LK));
% plot(TBL_M1.Depth_uM(IX_LK), Bfit, '-r')

% mean for Q_base and Q_post 1, 2 & 3
for ii = 1:length(IX_Lesion_LM1)
    TBL.Q_base_std(ii,:) = std(TBL.Q_Base{ii,1})/10;
    TBL.Q_post1_std(ii,:) = std(TBL.Q_Post1{ii,1})/10; 
    TBL.Q_post2_std(ii,:) = std(TBL.Q_Post2{ii,1})/10; 
    TBL.Q_post3_std(ii,:) = std(TBL.Q_Post3{ii,1})/10; 
end

% mann whitney of baseline and post 1, 2 & 3
for ii = 1:length(IX_CM1)
    TBL.Q_base_post1_p(ii) = ranksum(TBL.Q_Base{ii,1}, TBL.Q_Post1{ii,1});
    TBL.Q_base_post2_p(ii) = ranksum(TBL.Q_Base{ii,1}, TBL.Q_Post2{ii,1});
    TBL.Q_base_post3_p(ii) = ranksum(TBL.Q_Base{ii,1}, TBL.Q_Post3{ii,1});
end

% zscore and percentile of baseline and post 1, 2 & 3
for ii = 1:length(IX_CM1)
    TBL.Q_base_post1_p(ii) = ranksum(TBL.Q_Base{ii,1}, TBL.Q_Post1{ii,1});
    TBL.Q_base_post2_p(ii) = ranksum(TBL.Q_Base{ii,1}, TBL.Q_Post2{ii,1});
    TBL.Q_base_post3_p(ii) = ranksum(TBL.Q_Base{ii,1}, TBL.Q_Post3{ii,1});
end

x = categorical({'Post Sal' 'Post Ket Early' 'Post Ket Late'});
x = reordercats(x,{'Post Sal' 'Post Ket Early' 'Post Ket Late'});
figure;bar(x,[sum(TBL.Q_base_post1_p(IX_CM1)<.05) sum(TBL.Q_base_post2_p(IX_CM1)<.05) sum(TBL.Q_base_post3_p(IX_CM1)<.05)])
xlabel('Relative to baseline')
ylable('Neurons')
title('Mann Whitney U test <.05 total neurons = 208 Control M1')

% zscore of post 1, 2 & 3 and the baseline for each neuron
for ii = 1:length(IX_Lesion_LM1)
    TBL.Q_base_post1_z(ii) = (TBL.Frate_post1(ii)-TBL.Frate_base(ii))/TBL.Q_base_std(ii);
    TBL.Q_base_post2_z(ii) = (TBL.Frate_post2(ii)-TBL.Frate_base(ii))/TBL.Q_base_std(ii);
    TBL.Q_base_post3_z(ii) = (TBL.Frate_post3(ii)-TBL.Frate_base(ii))/TBL.Q_base_std(ii);
end
% changing the infs to nans in the z scores
find_inf_1(:,1) = find(TBL.Q_base_post1_z == inf);
TBL.Q_base_post1_z(find_inf_1)= nan;

find_inf_2(:,1) = find(TBL.Q_base_post2_z == inf);
TBL.Q_base_post2_z(find_inf_2)= nan;

find_inf_3(:,1) = find(TBL.Q_base_post3_z == inf);
TBL.Q_base_post3_z(find_inf_3)= nan;

% looking at which neurons increased/decreased post inj for
% neurons with zscores above/below .5/-.5
% Control M1
% x1 = categorical({'Increase FR' 'Decrease FR'});
% x1 = reordercats(x1,{'Increase FR' 'Decrease FR'});
% figure;bar(x1,[sum(IZ_inc_post2_CM1) sum(IZ_dec_post2_CM1)])
% xlabel('20 min Post ketamine')
% ylabel('Neurons')
% title('Zscores of FR (above/below .5/-.5) increase vs decrease Post ket Control M1')
% 
% for post sal and post 2
% Control M1
IZ_inc_postsal_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post1_z > .5;
IZ_dec_postsal_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post1_z < -.5;
IZ_nodiff_postsal_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;

IZ_inc_post2_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post2_z > .5;
IZ_dec_post2_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post2_z < -.5;
IZ_nodiff_post2_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;

IZ_inc_post3_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post3_z > .5;
IZ_dec_post3_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post3_z < -.5;
IZ_nodiff_post3_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;

x2 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x2 = reordercats(x2,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x2,[sum(IZ_inc_postsal_CM1) sum(IZ_dec_postsal_CM1) sum(IZ_nodiff_postsal_CM1); sum(IZ_inc_post2_CM1) sum(IZ_dec_post2_CM1) sum(IZ_nodiff_post2_CM1); sum(IZ_inc_post3_CM1) sum(IZ_dec_post3_CM1) sum(IZ_nodiff_post3_CM1)])
ylabel('Neurons')
title('Zscores of FR (above/below .5/-.5) increase vs decrease vs below cutoff Post inj Control M1')

y = [TBL.Q_base_post2_z(IZ_inc_post2_CM1);TBL.Q_base_post2_z(IZ_dec_post2_CM1)];
g1 = repmat({'Increased FR'},79,1);
g2 = repmat({'Decreased FR'},62,1);
g = [g1; g2];
figure;boxplot(y,g)
xlabel('20 min Post ketamine')
ylabel('zscore FR')
title('Zscores of FR (above/below .5/-.5) increase vs decrease Post ket Control M1')
% 6ODHA_LID M1
% x1 = categorical({'Increase FR' 'Decrease FR'});
% x1 = reordercats(x1,{'Increase FR' 'Decrease FR'});
% figure;bar(x1,[sum(IZ_inc_post2_LM1) sum(IZ_dec_post2_LM1)])
% xlabel('20 min Post ketamine')
% ylabel('Neurons')
% title('Zscores of FR (above/below .5/-.5) increase vs decrease Post ket LID M1')
%6ODHA_LID M1
IZ_inc_postsal_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post1_z > .5;
IZ_dec_postsal_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post1_z < -.5;
IZ_nodiff_postsal_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;

IZ_inc_post2_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post2_z > .5;
IZ_dec_post2_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post2_z < -.5;
IZ_nodiff_post2_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;

IZ_inc_post3_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post3_z > .5;
IZ_dec_post3_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post3_z < -.5;
IZ_nodiff_post3_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;

x3 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x3 = reordercats(x3,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x3,[sum(IZ_inc_postsal_LM1) sum(IZ_dec_postsal_LM1) sum(IZ_nodiff_postsal_LM1); sum(IZ_inc_post2_LM1) sum(IZ_dec_post2_LM1) sum(IZ_nodiff_post2_LM1); sum(IZ_inc_post3_LM1) sum(IZ_dec_post3_LM1) sum(IZ_nodiff_post3_LM1)])
ylabel('Neurons')
title('Zscores of FR (above/below .5/-.5) increase vs decrease vs below cutoff Post inj LID M1')


y = [TBL.Q_base_post2_z(IZ_inc_post2_LM1);TBL.Q_base_post2_z(IZ_dec_post2_LM1)];
g1 = repmat({'Increased FR'},26,1);
g2 = repmat({'Decreased FR'},18,1);
g = [g1; g2];
figure;boxplot(y,g)
xlabel('20 min Post ketamine')
ylabel('zscore FR')
title('Zscores of FR (above/below .5/-.5) increase vs decrease Post ket LID M1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non responder pool and responder pool for chi square analysis between conditions for increasers and decreasers and non responders
% Post and pre ketamine injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IR_resp =  TBL.Q_base_post1_z > .5 | TBL.Q_base_post1_z < -.5 | TBL.Q_base_post2_z > .5 | TBL.Q_base_post2_z < -.5 | TBL.Q_base_post3_z > .5 | TBL.Q_base_post3_z < -.5;
IR_resp_inc = TBL.Q_base_post1_z > .5 | TBL.Q_base_post2_z > .5 | TBL.Q_base_post3_z > .5;
IR_resp_dec = TBL.Q_base_post1_z < -.5 | TBL.Q_base_post2_z < -.5 | TBL.Q_base_post3_z < -.5;
IT_post1_noresp = TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;
IT_post2_noresp = TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;
IT_post3_noresp = TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;
IR_noresp = IT_post1_noresp | IT_post2_noresp | IT_post3_noresp;
% Control M1
% post saline increase decrease and no response of the responders pool
IR_inc_postsal_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post1_z > .5;
IR_dec_postsal_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post1_z < -.5;
IR_noresp_postsal_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;
% post ket increase decrease and no response of the responders pool
IR_inc_post2_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post2_z > .5;
IR_dec_post2_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post2_z < -.5;
IR_noresp_post2_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;
% post ket late increase decrease and no response of the responders pool
IR_inc_post3_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post3_z > .5;
IR_dec_post3_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post3_z < -.5;
IR_noresp_post3_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;

x4 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x4 = reordercats(x4,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x4,[sum(IR_inc_postsal_CM1) sum(IR_dec_postsal_CM1) sum(IR_noresp_postsal_CM1); sum(IR_inc_post2_CM1) sum(IR_dec_post2_CM1) sum(IR_noresp_post2_CM1); sum(IR_inc_post3_CM1) sum(IR_dec_post3_CM1) sum(IR_noresp_post3_CM1)], 'stacked')
ylabel('Neurons')
title('For responding neurons Zscores of FR (above/below .5/-.5) increase vs decrease vs below cutoff Post inj Control M1')

% chi quared binomial test for within conditions
p_sal_inc = myBinomTest(sum(IR_inc_postsal_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_sal_dec = myBinomTest(sum(IR_dec_postsal_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_sal_noresp = myBinomTest(sum(IR_noresp_postsal_CM1), sum(IX_CM1 & IR_resp), 1/3);

p_ket_inc = myBinomTest(sum(IR_inc_post2_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_ket_dec = myBinomTest(sum(IR_dec_post2_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_ket_noresp = myBinomTest(sum(IR_noresp_post2_CM1), sum(IX_CM1 & IR_resp), 1/3);

p_post3_inc = myBinomTest(sum(IR_inc_post3_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_post3_dec = myBinomTest(sum(IR_dec_post3_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_post3_noresp = myBinomTest(sum(IR_noresp_post3_CM1), sum(IX_CM1 & IR_resp), 1/3);

% binomial test for between conditions
p_cond_sal_inc = myBinomTest(sum(IR_inc_postsal_CM1), sum(IX_CM1 & IR_resp_inc), 1/3);
p_cond_ket_inc = myBinomTest(sum(IR_inc_post2_CM1), sum(IX_CM1 & IR_resp_inc), 1/3);
p_cond_post3_inc = myBinomTest(sum(IR_inc_post3_CM1), sum(IX_CM1 & IR_resp_inc), 1/3);

p_cond_sal_dec = myBinomTest(sum(IR_dec_postsal_CM1), sum(IX_CM1 & IR_resp_dec), 1/3);
p_cond_ket_dec = myBinomTest(sum(IR_dec_post2_CM1), sum(IX_CM1 & IR_resp_dec), 1/3);
p_cond_post3_dec = myBinomTest(sum(IR_dec_post3_CM1), sum(IX_CM1 & IR_resp_dec), 1/3);

p_cond_sal_noresp = myBinomTest(sum(IR_noresp_postsal_CM1), sum(IX_CM1 & IR_noresp), 1/3);
p_cond_ket_noresp = myBinomTest(sum(IR_noresp_post2_CM1), sum(IX_CM1 & IR_noresp), 1/3);
p_cond_post3_noresp = myBinomTest(sum(IR_noresp_post3_CM1), sum(IX_CM1 & IR_noresp), 1/3);

% Control Straitum
IR_inc_postsal_CS = IX_CS & IR_resp & TBL.Q_base_post1_z > .5;
IR_dec_postsal_CS = IX_CS & IR_resp & TBL.Q_base_post1_z < -.5;
IR_noresp_postsal_CS = IX_CS & IR_resp & TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;
% post ket increase decrease and no response of the responders pool
IR_inc_post2_CS = IX_CS & IR_resp & TBL.Q_base_post2_z > .5;
IR_dec_post2_CS = IX_CS & IR_resp & TBL.Q_base_post2_z < -.5;
IR_noresp_post2_CS = IX_CS & IR_resp & TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;
% post ket late increase decrease and no response of the responders pool
IR_inc_post3_CS = IX_CS & IR_resp & TBL.Q_base_post3_z > .5;
IR_dec_post3_CS = IX_CS & IR_resp & TBL.Q_base_post3_z < -.5;
IR_noresp_post3_CS = IX_CS & IR_resp & TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;

x4 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x4 = reordercats(x4,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x4,[sum(IR_inc_postsal_CS) sum(IR_dec_postsal_CS) sum(IR_noresp_postsal_CS); sum(IR_inc_post2_CS) ...
    sum(IR_dec_post2_CS) sum(IR_noresp_post2_CS); sum(IR_inc_post3_CS) sum(IR_dec_post3_CS) sum(IR_noresp_post3_CS)], 'stacked')
ylabel('Neurons')
title('For responding neurons Zscores of FR (above/below .5/-.5) increase vs decrease vs below cutoff Post inj Control Striatum')


% 6ODHA_LID M1
%%%%%%% lesioned hemisphere %%%%%%%%%%%%%
% post saline increase decrease and no response of the responders pool
IR_inc_postsal_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_z > .5;
IR_dec_postsal_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_z < -.5;
IR_noresp_postsal_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;
% post ket increase decrease and no response of the responders pool
IR_inc_post2_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_z > .5;
IR_dec_post2_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_z < -.5;
IR_noresp_post2_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;
% post ket late increase decrease and no response of the responders pool
IR_inc_post3_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_z > .5;
IR_dec_post3_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_z < -.5;
IR_noresp_post3_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;

x5 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x5 = reordercats(x5,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x5,[sum(IR_inc_postsal_LM1) sum(IR_dec_postsal_LM1) sum(IR_noresp_postsal_LM1); sum(IR_inc_post2_LM1) ... 
    sum(IR_dec_post2_LM1) sum(IR_noresp_post2_LM1); sum(IR_inc_post3_LM1) sum(IR_dec_post3_LM1) sum(IR_noresp_post3_LM1)], 'stacked')
ylabel('Neurons')
title('For responding neurons Zscores of FR (above/below .5/-.5) increase vs decrease vs below cutoff Post inj LID M1 Lesioned hemsiphere')

% chi quared binomial test for within conditions
p_sal_inc_LM1 = myBinomTest(sum(IR_inc_postsal_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_sal_dec_LM1 = myBinomTest(sum(IR_dec_postsal_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_sal_noresp_LM1 = myBinomTest(sum(IR_noresp_postsal_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

p_ket_inc_LM1 = myBinomTest(sum(IR_inc_post2_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_ket_dec_LM1 = myBinomTest(sum(IR_dec_post2_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_ket_noresp_LM1 = myBinomTest(sum(IR_noresp_post2_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

p_post3_inc_LM1 = myBinomTest(sum(IR_inc_post3_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_post3_dec_LM1 = myBinomTest(sum(IR_dec_post3_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_post3_noresp_LM1 = myBinomTest(sum(IR_noresp_post3_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

% binomial test for between conditions
p_cond_sal_inc_LM1 = myBinomTest(sum(IR_inc_postsal_LM1), sum(IX_Lesion_LM1 & IR_resp_inc), 1/3);
p_cond_ket_inc_LM1 = myBinomTest(sum(IR_inc_post2_LM1), sum(IX_Lesion_LM1 & IR_resp_inc), 1/3);
p_cond_post3_inc_LM1 = myBinomTest(sum(IR_inc_post3_LM1), sum(IX_Lesion_LM1 & IR_resp_inc), 1/3);

p_cond_sal_dec_LM1 = myBinomTest(sum(IR_dec_postsal_LM1), sum(IX_Lesion_LM1 & IR_resp_dec), 1/3);
p_cond_ket_dec_LM1 = myBinomTest(sum(IR_dec_post2_LM1), sum(IX_Lesion_LM1 & IR_resp_dec), 1/3);
p_cond_post3_dec_LM1 = myBinomTest(sum(IR_dec_post3_LM1), sum(IX_Lesion_LM1 & IR_resp_dec), 1/3);

p_cond_sal_noresp_LM1 = myBinomTest(sum(IR_noresp_postsal_LM1), sum(IX_Lesion_LM1 & IR_noresp), 1/3);
p_cond_ket_noresp_LM1 = myBinomTest(sum(IR_noresp_post2_LM1), sum(IX_Lesion_LM1 & IR_noresp), 1/3);
p_cond_post3_noresp_LM1 = myBinomTest(sum(IR_noresp_post3_LM1), sum(IX_Lesion_LM1 & IR_noresp), 1/3);

%%%%%% Unlesioend hemisphere %%%%%%%%%%%
IR_inc_postsal_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_z > .5;
IR_dec_postsal_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_z < -.5;
IR_noresp_postsal_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;
% post ket increase decrease and no response of the responders pool
IR_inc_post2_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_z > .5;
IR_dec_post2_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_z < -.5;
IR_noresp_post2_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;
% post ket late increase decrease and no response of the responders pool
IR_inc_post3_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_z > .5;
IR_dec_post3_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_z < -.5;
IR_noresp_post3_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;

x5 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x5 = reordercats(x5,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x5,[sum(IR_inc_postsal_LM1) sum(IR_dec_postsal_LM1) sum(IR_noresp_postsal_LM1); sum(IR_inc_post2_LM1) ...
    sum(IR_dec_post2_LM1) sum(IR_noresp_post2_LM1); sum(IR_inc_post3_LM1) sum(IR_dec_post3_LM1) sum(IR_noresp_post3_LM1)], 'stacked')
ylabel('Neurons')
title('For responding neurons Zscores of FR (above/below .5/-.5) increase vs decrease vs below cutoff Post inj LID M1 Unlesioned')

% chi quared binomial test for within conditions
p_sal_inc_LM1 = myBinomTest(sum(IR_inc_postsal_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_sal_dec_LM1 = myBinomTest(sum(IR_dec_postsal_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_sal_noresp_LM1 = myBinomTest(sum(IR_noresp_postsal_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

p_ket_inc_LM1 = myBinomTest(sum(IR_inc_post2_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_ket_dec_LM1 = myBinomTest(sum(IR_dec_post2_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_ket_noresp_LM1 = myBinomTest(sum(IR_noresp_post2_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

p_post3_inc_LM1 = myBinomTest(sum(IR_inc_post3_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_post3_dec_LM1 = myBinomTest(sum(IR_dec_post3_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_post3_noresp_LM1 = myBinomTest(sum(IR_noresp_post3_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

% binomial test for between conditions
p_cond_sal_inc_LM1 = myBinomTest(sum(IR_inc_postsal_LM1), sum(IX_UnLesion_LM1 & IR_resp_inc), 1/3);
p_cond_ket_inc_LM1 = myBinomTest(sum(IR_inc_post2_LM1), sum(IX_UnLesion_LM1 & IR_resp_inc), 1/3);
p_cond_post3_inc_LM1 = myBinomTest(sum(IR_inc_post3_LM1), sum(IX_UnLesion_LM1 & IR_resp_inc), 1/3);

p_cond_sal_dec_LM1 = myBinomTest(sum(IR_dec_postsal_LM1), sum(IX_UnLesion_LM1 & IR_resp_dec), 1/3);
p_cond_ket_dec_LM1 = myBinomTest(sum(IR_dec_post2_LM1), sum(IX_UnLesion_LM1 & IR_resp_dec), 1/3);
p_cond_post3_dec_LM1 = myBinomTest(sum(IR_dec_post3_LM1), sum(IX_UnLesion_LM1 & IR_resp_dec), 1/3);

p_cond_sal_noresp_LM1 = myBinomTest(sum(IR_noresp_postsal_LM1), sum(IX_UnLesion_LM1 & IR_noresp), 1/3);
p_cond_ket_noresp_LM1 = myBinomTest(sum(IR_noresp_post2_LM1), sum(IX_UnLesion_LM1 & IR_noresp), 1/3);
p_cond_post3_noresp_LM1 = myBinomTest(sum(IR_noresp_post3_LM1), sum(IX_UnLesion_LM1 & IR_noresp), 1/3);

%% Looking at LDOPA pre and post change in z score firing rates
% 6ODHA_LID M1
%%%%%%%%%%%%%% Lesioned Hemsiphere %%%%%%%%%%%%%%%%%%%%%%
% post saline increase decrease and no response of the responders pool
IR_inc_postLDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_z > .5;
IR_dec_postLDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_z < -.5;
IR_noresp_postLDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;
% post ket increase decrease and no response of the responders pool
IR_inc_postLDO_80_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_z > .5;
IR_dec_postLDO_80_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_z < -.5;
IR_noresp_postLDO_80_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;
% post ket late increase decrease and no response of the responders pool
IR_inc_post3_LDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_z > .5;
IR_dec_post3_LDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_z < -.5;
IR_noresp_post3_LDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;

x6 = categorical({'Post LDOPA' 'Peak 80' 'Post 60 min'});
x6 = reordercats(x6,{'Post LDOPA' 'Peak 80' 'Post 60 min'});
figure;bar(x6,[sum(IR_inc_postLDO_LM1) sum(IR_dec_postLDO_LM1) sum(IR_noresp_postLDO_LM1); sum(IR_inc_postLDO_80_LM1) ...
    sum(IR_dec_postLDO_80_LM1) sum(IR_noresp_postLDO_80_LM1); sum(IR_inc_post3_LDO_LM1) sum(IR_dec_post3_LDO_LM1) sum(IR_noresp_post3_LDO_LM1)], 'stacked')
ylabel('Neurons')
title('For responding neurons Zscores of FR (above/below .5/-.5) increase vs decrease vs below cutoff Post LDOPA inj LID M1 lesioned hemisphere')

% chi quared binomial test for within conditions
p_LDO_inc_LM1 = myBinomTest(sum(IR_inc_postLDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_dec_LM1 = myBinomTest(sum(IR_dec_postLDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_noresp_LM1 = myBinomTest(sum(IR_noresp_postLDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

p_LDO_80_inc_LM1 = myBinomTest(sum(IR_inc_postLDO_80_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_80_dec_LM1 = myBinomTest(sum(IR_dec_postLDO_80_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_80_noresp_LM1 = myBinomTest(sum(IR_noresp_postLDO_80_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

p_post3_LDO_inc_LM1 = myBinomTest(sum(IR_inc_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_post3_LDO_dec_LM1 = myBinomTest(sum(IR_dec_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_post3_LDO_noresp_LM1 = myBinomTest(sum(IR_noresp_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

% binomial test for between conditions
p_cond_LDO_inc_LM1 = myBinomTest(sum(IR_inc_postLDO_LM1), sum(IX_Lesion_LM1 & IR_resp_inc), 1/3);
p_cond_LDO_80_inc_LM1 = myBinomTest(sum(IR_inc_postLDO_80_LM1), sum(IX_Lesion_LM1 & IR_resp_inc), 1/3);
p_cond_post3_LDO_inc_LM1 = myBinomTest(sum(IR_inc_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_resp_inc), 1/3);

p_cond_LDO_dec_LM1 = myBinomTest(sum(IR_dec_postLDO_LM1), sum(IX_Lesion_LM1 & IR_resp_dec), 1/3);
p_cond_LDO_80_dec_LM1 = myBinomTest(sum(IR_dec_postLDO_80_LM1), sum(IX_Lesion_LM1 & IR_resp_dec), 1/3);
p_cond_post3_LDO_dec_LM1 = myBinomTest(sum(IR_dec_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_resp_dec), 1/3);

p_cond_LDO_noresp_LM1 = myBinomTest(sum(IR_noresp_postLDO_LM1), sum(IX_Lesion_LM1 & IR_noresp), 1/3);
p_cond_LDO_80_noresp_LM1 = myBinomTest(sum(IR_noresp_postLDO_80_LM1), sum(IX_Lesion_LM1 & IR_noresp), 1/3);
p_cond_post3_LDO_noresp_LM1 = myBinomTest(sum(IR_noresp_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_noresp), 1/3);

%%%%%%%%%%%%%%%% Unlesioned hemisphere %%%%%%%%%%%%%%%
IR_inc_postLDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_z > .5;
IR_dec_postLDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_z < -.5;
IR_noresp_postLDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_z < .5 & TBL.Q_base_post1_z > -.5;
% post ket increase decrease and no response of the responders pool
IR_inc_postLDO_80_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_z > .5;
IR_dec_postLDO_80_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_z < -.5;
IR_noresp_postLDO_80_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_z < .5 & TBL.Q_base_post2_z > -.5;
% post ket late increase decrease and no response of the responders pool
IR_inc_post3_LDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_z > .5;
IR_dec_post3_LDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_z < -.5;
IR_noresp_post3_LDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_z < .5 & TBL.Q_base_post3_z > -.5;

x6 = categorical({'Post LDOPA' 'Peak 80' 'Post 60 min'});
x6 = reordercats(x6,{'Post LDOPA' 'Peak 80' 'Post 60 min'});
figure;bar(x6,[sum(IR_inc_postLDO_LM1) sum(IR_dec_postLDO_LM1) sum(IR_noresp_postLDO_LM1); sum(IR_inc_postLDO_80_LM1) sum(IR_dec_postLDO_80_LM1) ...
    sum(IR_noresp_postLDO_80_LM1); sum(IR_inc_post3_LDO_LM1) sum(IR_dec_post3_LDO_LM1) sum(IR_noresp_post3_LDO_LM1)], 'stacked')
ylabel('Neurons')
title('For responding neurons Zscores of FR (above/below .5/-.5) increase vs decrease vs below cutoff Post LDOPA inj LID M1 Unlesioned hemisphere')

% chi quared binomial test for within conditions
p_LDO_inc_LM1 = myBinomTest(sum(IR_inc_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_dec_LM1 = myBinomTest(sum(IR_dec_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_noresp_LM1 = myBinomTest(sum(IR_noresp_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

p_LDO_80_inc_LM1 = myBinomTest(sum(IR_inc_postLDO_80_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_80_dec_LM1 = myBinomTest(sum(IR_dec_postLDO_80_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_80_noresp_LM1 = myBinomTest(sum(IR_noresp_postLDO_80_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

p_post3_LDO_inc_LM1 = myBinomTest(sum(IR_inc_post3_LDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_post3_LDO_dec_LM1 = myBinomTest(sum(IR_dec_post3_LDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_post3_LDO_noresp_LM1 = myBinomTest(sum(IR_noresp_post3_LDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

% binomial test for between conditions
p_cond_LDO_inc_LM1 = myBinomTest(sum(IR_inc_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_resp_inc), 1/3);
p_cond_LDO_80_inc_LM1 = myBinomTest(sum(IR_inc_postLDO_80_LM1), sum(IX_UnLesion_LM1 & IR_resp_inc), 1/3);
p_cond_post3_LDO_inc_LM1 = myBinomTest(sum(IR_inc_post3_LDO_LM1), sum(IX_UnLesion_LM1 & IR_resp_inc), 1/3);

p_cond_LDO_dec_LM1 = myBinomTest(sum(IR_dec_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_resp_dec), 1/3);
p_cond_LDO_80_dec_LM1 = myBinomTest(sum(IR_dec_postLDO_80_LM1), sum(IX_UnLesion_LM1 & IR_resp_dec), 1/3);
p_cond_post3_LDO_dec_LM1 = myBinomTest(sum(IR_dec_post3_LDO_LM1), sum(IX_UnLesion_LM1 & IR_resp_dec), 1/3);

p_cond_LDO_noresp_LM1 = myBinomTest(sum(IR_noresp_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_noresp), 1/3);
p_cond_LDO_80_noresp_LM1 = myBinomTest(sum(IR_noresp_postLDO_80_LM1), sum(IX_UnLesion_LM1 & IR_noresp), 1/3);
p_cond_post3_LDO_noresp_LM1 = myBinomTest(sum(IR_noresp_post3_LDO_LM1), sum(IX_UnLesion_LM1 & IR_noresp), 1/3);

%% difference between post 1, 2 & 3 with baseline mFR

for ii = 1:length(IX_Lesion_LM1)
    TBL.Q_base_post1_diff_mFR(ii) = (TBL.Frate_post1(ii)-TBL.Frate_base(ii));
    TBL.Q_base_post2_diff_mFR(ii) = (TBL.Frate_post2(ii)-TBL.Frate_base(ii));
    TBL.Q_base_post3_diff_mFR(ii) = (TBL.Frate_post3(ii)-TBL.Frate_base(ii));
end
% Getting the responders 
IR_resp = TBL.Q_base_post1_diff_mFR > .5| TBL.Q_base_post1_diff_mFR < -.5 | TBL.Q_base_post2_diff_mFR > .5 | TBL.Q_base_post2_diff_mFR < -.5 ...
    | TBL.Q_base_post3_diff_mFR > .5 | TBL.Q_base_post3_diff_mFR < -.5;

% Looking at the FR diff to see if it increased, decreased or did not
% change
%%% Pre and Post KETAMINE %%%
% Control M1
ID_inc_post2_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post2_diff_mFR > .5;
ID_dec_post2_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post2_diff_mFR < -.5;
ID_nodiff_post2_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post2_diff_mFR > -.5 & TBL.Q_base_post2_diff_mFR < .5;

ID_inc_postsal_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post1_diff_mFR > .5;
ID_dec_postsal_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post1_diff_mFR < -.5;
ID_nodiff_postsal_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post1_diff_mFR > -.5 & TBL.Q_base_post1_diff_mFR < .5;

ID_inc_post3_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post3_diff_mFR > .5;
ID_dec_post3_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post3_diff_mFR < -.5;
ID_nodiff_post3_CM1 = IX_CM1 & IR_resp & TBL.Q_base_post3_diff_mFR > -.5 & TBL.Q_base_post3_diff_mFR < .5;

x3 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x3 = reordercats(x3,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x3,[sum(ID_inc_postsal_CM1) sum(ID_dec_postsal_CM1) sum(ID_nodiff_postsal_CM1); sum(ID_inc_post2_CM1) sum(ID_dec_post2_CM1) sum(ID_nodiff_post2_CM1); sum(ID_inc_post3_CM1) sum(ID_dec_post3_CM1) sum(ID_nodiff_post3_CM1)], 'stacked')
ylabel('Neurons')
title('Difference in mean FR increase vs decrease vs no change (cutoff .5) Post inj Control M1')


% chi quared binomial test for within conditions
p_sal_inc = myBinomTest(sum(ID_inc_postsal_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_sal_dec = myBinomTest(sum(ID_dec_postsal_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_sal_noresp = myBinomTest(sum(ID_nodiff_postsal_CM1), sum(IX_CM1 & IR_resp), 1/3);

p_ket_inc = myBinomTest(sum(ID_inc_post2_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_ket_dec = myBinomTest(sum(ID_dec_post2_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_ket_noresp = myBinomTest(sum(ID_nodiff_post2_CM1), sum(IX_CM1 & IR_resp), 1/3);

p_post3_inc = myBinomTest(sum(ID_inc_post3_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_post3_dec = myBinomTest(sum(ID_dec_post3_CM1), sum(IX_CM1 & IR_resp), 1/3);
p_post3_noresp = myBinomTest(sum(ID_nodiff_post3_CM1), sum(IX_CM1 & IR_resp), 1/3);

% binomial test for between conditions
p_cond_sal_inc = myBinomTest(sum(ID_inc_postsal_CM1), sum(IX_CM1 & IR_resp_inc), 1/3);
p_cond_ket_inc = myBinomTest(sum(ID_inc_post2_CM1), sum(IX_CM1 & IR_resp_inc), 1/3);
p_cond_post3_inc = myBinomTest(sum(ID_inc_post3_CM1), sum(IX_CM1 & IR_resp_inc), 1/3);

p_cond_sal_dec = myBinomTest(sum(ID_dec_postsal_CM1), sum(IX_CM1 & IR_resp_dec), 1/3);
p_cond_ket_dec = myBinomTest(sum(ID_dec_post2_CM1), sum(IX_CM1 & IR_resp_dec), 1/3);
p_cond_post3_dec = myBinomTest(sum(ID_dec_post3_CM1), sum(IX_CM1 & IR_resp_dec), 1/3);

p_cond_sal_noresp = myBinomTest(sum(ID_nodiff_postsal_CM1), sum(IX_CM1 & IR_noresp), 1/3);
p_cond_ket_noresp = myBinomTest(sum(ID_nodiff_post2_CM1), sum(IX_CM1 & IR_noresp), 1/3);
p_cond_post3_noresp = myBinomTest(sum(IR_noresp_post3_CM1), sum(IX_CM1 & IR_noresp), 1/3);

% Control Striatum
ID_inc_post2_CS = IX_CS & IR_resp & TBL.Q_base_post2_diff_mFR > .5;
ID_dec_post2_CS = IX_CS & IR_resp & TBL.Q_base_post2_diff_mFR < -.5;
ID_nodiff_post2_CS = IX_CS & IR_resp & TBL.Q_base_post2_diff_mFR > -.5 & TBL.Q_base_post2_diff_mFR < .5;

ID_inc_postsal_CS = IX_CS & IR_resp & TBL.Q_base_post1_diff_mFR > .5;
ID_dec_postsal_CS = IX_CS & IR_resp & TBL.Q_base_post1_diff_mFR < -.5;
ID_nodiff_postsal_CS = IX_CS & IR_resp & TBL.Q_base_post1_diff_mFR > -.5 & TBL.Q_base_post1_diff_mFR < .5;

ID_inc_post3_CS = IX_CS & IR_resp & TBL.Q_base_post3_diff_mFR > .5;
ID_dec_post3_CS = IX_CS & IR_resp & TBL.Q_base_post3_diff_mFR < -.5;
ID_nodiff_post3_CS = IX_CS & IR_resp & TBL.Q_base_post3_diff_mFR > -.5 & TBL.Q_base_post3_diff_mFR < .5;

x3 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x3 = reordercats(x3,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x3,[sum(ID_inc_postsal_CS) sum(ID_dec_postsal_CS) sum(ID_nodiff_postsal_CS); sum(ID_inc_post2_CS) ...
    sum(ID_dec_post2_CS) sum(ID_nodiff_post2_CS); sum(ID_inc_post3_CS) sum(ID_dec_post3_CS) sum(ID_nodiff_post3_CS)], 'stacked')
ylabel('Neurons')
title('Difference in mean FR increase vs decrease vs no change (cutoff .5) Post inj Control Striatum')

% chi quared binomial test for within conditions
p_sal_inc = myBinomTest(sum(ID_inc_postsal_CS), sum(IX_CS & IR_resp), 1/3);
p_sal_dec = myBinomTest(sum(ID_dec_postsal_CS), sum(IX_CS & IR_resp), 1/3);
p_sal_noresp = myBinomTest(sum(ID_nodiff_postsal_CS), sum(IX_CS & IR_resp), 1/3);

p_ket_inc = myBinomTest(sum(ID_inc_post2_CS), sum(IX_CS & IR_resp), 1/3);
p_ket_dec = myBinomTest(sum(ID_dec_post2_CS), sum(IX_CS & IR_resp), 1/3);
p_ket_noresp = myBinomTest(sum(ID_nodiff_post2_CS), sum(IX_CS & IR_resp), 1/3);

p_post3_inc = myBinomTest(sum(ID_inc_post3_CS), sum(IX_CS & IR_resp), 1/3);
p_post3_dec = myBinomTest(sum(ID_dec_post3_CS), sum(IX_CS & IR_resp), 1/3);
p_post3_noresp = myBinomTest(sum(ID_nodiff_post3_CS), sum(IX_CS & IR_resp), 1/3);



% 6ODHA_LID
%%%%%%%%%% lesioned hemisphere %%%%%%%
ID_inc_post2_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR > .5;
ID_dec_post2_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR < -.5;
ID_nodiff_post2_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR > -.5 & TBL.Q_base_post2_diff_mFR < .5;

ID_inc_postsal_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR > .5;
ID_dec_postsal_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR < -.5;
ID_nodiff_postsal_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR > -.5 & TBL.Q_base_post1_diff_mFR < .5;

ID_inc_post3_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR > .5;
ID_dec_post3_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR < -.5;
ID_nodiff_post3_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR > -.5 & TBL.Q_base_post3_diff_mFR < .5;

x3 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x3 = reordercats(x3,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x3,[sum(ID_inc_postsal_LM1) sum(ID_dec_postsal_LM1) sum(ID_nodiff_postsal_LM1); sum(ID_inc_post2_LM1) sum(ID_dec_post2_LM1) sum(ID_nodiff_post2_LM1); sum(ID_inc_post3_LM1) sum(ID_dec_post3_LM1) sum(ID_nodiff_post3_LM1)], 'stacked')
set(gca,'ylim',[0 35])
ylabel('Neurons')
title('Difference in mean FR increase vs decrease vs no change Post inj Lesioned LID M1')

% chi quared binomial test for within conditions
p_sal_inc_LM1 = myBinomTest(sum(ID_inc_postsal_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_sal_dec_LM1 = myBinomTest(sum(ID_dec_postsal_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_sal_noresp_LM1 = myBinomTest(sum(ID_nodiff_postsal_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

p_ket_inc_LM1 = myBinomTest(sum(ID_inc_post2_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_ket_dec_LM1 = myBinomTest(sum(ID_dec_post2_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_ket_noresp_LM1 = myBinomTest(sum(ID_nodiff_post2_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

p_post3_inc_LM1 = myBinomTest(sum(ID_inc_post3_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_post3_dec_LM1 = myBinomTest(sum(ID_dec_post3_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_post3_noresp_LM1 = myBinomTest(sum(ID_nodiff_post3_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

%%%%%%%%%% Unlesioned hem %%%%%%%%%%%%%%
ID_inc_post2_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR > .5;
ID_dec_post2_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR < -.5;
ID_nodiff_post2_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR > -.5 & TBL.Q_base_post2_diff_mFR < .5;

ID_inc_postsal_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR > .5;
ID_dec_postsal_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR < -.5;
ID_nodiff_postsal_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR > -.5 & TBL.Q_base_post1_diff_mFR < .5;

ID_inc_post3_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR > .5;
ID_dec_post3_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR < -.5;
ID_nodiff_post3_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR > -.5 & TBL.Q_base_post3_diff_mFR < .5;

x3 = categorical({'Post sal' 'Post Ket' 'Post Ket late'});
x3 = reordercats(x3,{'Post sal' 'Post Ket' 'Post Ket late'});
figure;bar(x3,[sum(ID_inc_postsal_LM1) sum(ID_dec_postsal_LM1) sum(ID_nodiff_postsal_LM1); sum(ID_inc_post2_LM1) sum(ID_dec_post2_LM1) sum(ID_nodiff_post2_LM1); sum(ID_inc_post3_LM1) sum(ID_dec_post3_LM1) sum(ID_nodiff_post3_LM1)],'stacked')
ylabel('Neurons')
title('Difference in mean FR increase vs decrease vs no change Post inj Unlesioned LID M1')

% chi quared binomial test for within conditions
%% Looking at FR diff pre and post LDOPA %%%
% 6ODHA_LID
%%%%%%%% Lesioned hemisphere %%%%%%%%%%
ID_inc_LDO_80_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR > .5;
ID_dec_LDO_80_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR < -.5;
ID_nodiff_LDO_80_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR > -.5 & TBL.Q_base_post2_diff_mFR < .5;

ID_inc_postLDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR > .5;
ID_dec_postLDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR < -.5;
ID_nodiff_postLDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR > -.5 & TBL.Q_base_post1_diff_mFR < .5;

ID_inc_post3_LDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR > .5;
ID_dec_post3_LDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR < -.5;
ID_nodiff_post3_LDO_LM1 = IX_Lesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR > -.5 & TBL.Q_base_post3_diff_mFR < .5;

x7 = categorical({'Post LDOPA' 'Peak 80' 'Post 60 min'});
x7 = reordercats(x7,{'Post LDOPA' 'Peak 80' 'Post 60 min'});
figure;bar(x7,[sum(ID_inc_postLDO_LM1) sum(ID_dec_postLDO_LM1) sum(ID_nodiff_postLDO_LM1); sum(ID_inc_LDO_80_LM1) ...
    sum(ID_dec_LDO_80_LM1) sum(ID_nodiff_LDO_80_LM1); sum(ID_inc_post3_LDO_LM1) sum(ID_dec_post3_LDO_LM1) sum(ID_nodiff_post3_LDO_LM1)], 'stacked')
ylabel('Neurons')
title('Difference in mean FR increase vs decrease vs no change (.5 cutoff) Post LDOPA inj LID M1 lesion hemisphere')
% testing
p_LDO_inc_LM1(1) = myBinomTest(sum(ID_inc_postLDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_dec_LM1(1) = myBinomTest(sum(ID_dec_postLDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_noresp_LM1(1) = myBinomTest(sum(ID_nodiff_postLDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

p_LDO_inc_LM1(2) = myBinomTest(sum(ID_inc_LDO_80_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_dec_LM1(2) = myBinomTest(sum(ID_dec_LDO_80_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_noresp_LM1(2) = myBinomTest(sum(ID_nodiff_LDO_80_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

p_LDO_inc_LM1(3) = myBinomTest(sum(ID_inc_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_dec_LM1(3) = myBinomTest(sum(ID_dec_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);
p_LDO_noresp_LM1(3) = myBinomTest(sum(ID_nodiff_post3_LDO_LM1), sum(IX_Lesion_LM1 & IR_resp), 1/3);

%%%%%%%%%%%%% Unlesioned hemisphere %%%%%%%%%%%%%
ID_inc_LDO_80_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR > .5;
ID_dec_LDO_80_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR < -.5;
ID_nodiff_LDO_80_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post2_diff_mFR > -.5 & TBL.Q_base_post2_diff_mFR < .5;

ID_inc_postLDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR > .5;
ID_dec_postLDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR < -.5;
ID_nodiff_postLDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post1_diff_mFR > -.5 & TBL.Q_base_post1_diff_mFR < .5;

ID_inc_post3_LDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR > .5;
ID_dec_post3_LDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR < -.5;
ID_nodiff_post3_LDO_LM1 = IX_UnLesion_LM1 & IR_resp & TBL.Q_base_post3_diff_mFR > -.5 & TBL.Q_base_post3_diff_mFR < .5;

x8 = categorical({'Post LDOPA' 'Peak 80' 'Post 60 min'});
x8 = reordercats(x8,{'Post LDOPA' 'Peak 80' 'Post 60 min'});
figure;bar(x8,[sum(ID_inc_postLDO_LM1) sum(ID_dec_postLDO_LM1) sum(ID_nodiff_postLDO_LM1); sum(ID_inc_LDO_80_LM1) sum(ID_dec_LDO_80_LM1) ...
    sum(ID_nodiff_LDO_80_LM1); sum(ID_inc_post3_LDO_LM1) sum(ID_dec_post3_LDO_LM1) sum(ID_nodiff_post3_LDO_LM1)], 'stacked')
ylabel('Neurons')
title('Difference in mean FR increase vs decrease vs no change (.5 cutoff) Post LDOPA inj LID M1 Unlesion hemisphere')

% Chi squard test
p_LDO_inc_LM1(1) = myBinomTest(sum(ID_inc_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_dec_LM1(1) = myBinomTest(sum(ID_dec_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_noresp_LM1(1) = myBinomTest(sum(ID_nodiff_postLDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

p_LDO_inc_LM1(2) = myBinomTest(sum(ID_inc_LDO_80_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_dec_LM1(2) = myBinomTest(sum(ID_dec_LDO_80_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_noresp_LM1(2) = myBinomTest(sum(ID_nodiff_LDO_80_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

p_LDO_inc_LM1(3) = myBinomTest(sum(ID_inc_post3_LDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_dec_LM1(3) = myBinomTest(), sum(IX_UnLesion_LM1 & IR_resp), 1/3);
p_LDO_noresp_LM1(3) = myBinomTest(sum(ID_nodiff_post3_LDO_LM1), sum(IX_UnLesion_LM1 & IR_resp), 1/3);

%% paired test on the diff
[h,p_post1_post2_CM1] = ttest(TBL.Q_base_post1_diff_mFR(IX_CM1),TBL.Q_base_post2_diff_mFR(IX_CM1));
[h,p_post1_post3_CM1] = ttest(TBL.Q_base_post1_diff_mFR(IX_CM1),TBL.Q_base_post3_diff_mFR(IX_CM1));
[h,p_post2_post3_CM1] = ttest(TBL.Q_base_post2_diff_mFR(IX_CM1),TBL.Q_base_post3_diff_mFR(IX_CM1));

% boxplot of diff
figure
boxplot([TBL.Q_base_post1_diff_mFR(IX_CM1) TBL.Q_base_post2_diff_mFR(IX_CM1) TBL.Q_base_post3_diff_mFR(IX_CM1)])
yt = get(gca, 'YTick');
axis([xlim    -200  ceil(max(yt)*1.5)])
xt = get(gca, 'XTick');
title('Diff from baseline; paired ttest')
hold on
if p_post1_post2_CM1<.05
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
if p_post2_post3_CM1<.05
    plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k')
end
if p_post1_post3_CM1<.05
    plot(xt([1 3]), [1 1]*max(yt)*1.3, '-k',  mean(xt([1 3])), max(yt)*1.35, '*k')
end
hold off

% paired test on the diff for Control Striatum
[h,p_post1_post2_CS] = ttest(TBL.Q_base_post1_diff_mFR(IX_CS),TBL.Q_base_post2_diff_mFR(IX_CS));
[h,p_post1_post3_CS] = ttest(TBL.Q_base_post1_diff_mFR(IX_CS),TBL.Q_base_post3_diff_mFR(IX_CS));
[h,p_post2_post3_CS] = ttest(TBL.Q_base_post2_diff_mFR(IX_CS),TBL.Q_base_post3_diff_mFR(IX_CS));

% boxplot of diff
figure
boxplot([TBL.Q_base_post1_diff_mFR(IX_CS) TBL.Q_base_post2_diff_mFR(IX_CS) TBL.Q_base_post3_diff_mFR(IX_CS)])
yt = get(gca, 'YTick');
axis([xlim    -20  ceil(max(yt)*1.5)])
xt = get(gca, 'XTick');
title('Diff from baseline Control Striatum; paired ttest')
hold on
if p_post1_post2_CS<.05
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
if p_post2_post3_CS<.05
    plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k')
end
if p_post1_post3_CS<.05
    plot(xt([1 3]), [1 1]*max(yt)*1.3, '-k',  mean(xt([1 3])), max(yt)*1.35, '*k')
end
hold off

%%%%%%%%%%% Lesioned hemisphere %%%%%%%%%%%%
% paired test on the diff for LID M1
[h,p_post1_post2_LM1] = ttest(TBL.Q_base_post1_diff_mFR(IX_Lesion_LM1),TBL.Q_base_post2_diff_mFR(IX_Lesion_LM1));
[h,p_post1_post3_LM1] = ttest(TBL.Q_base_post1_diff_mFR(IX_Lesion_LM1),TBL.Q_base_post3_diff_mFR(IX_Lesion_LM1));
[h,p_post2_post3_LM1] = ttest(TBL.Q_base_post2_diff_mFR(IX_Lesion_LM1),TBL.Q_base_post3_diff_mFR(IX_Lesion_LM1));

% boxplot of diff
figure
violin([TBL.Q_base_post1_diff_mFR(IX_Lesion_LM1) TBL.Q_base_post2_diff_mFR(IX_Lesion_LM1) TBL.Q_base_post3_diff_mFR(IX_Lesion_LM1)], ...
    'xlabel',{'Post LDOPA','Peak 80','Post 60 min'})
yt = get(gca, 'YTick');
axis([xlim    -100  ceil(max(yt)*1.5)])
xt = get(gca, 'XTick');
title('Diff from baseline LID M1; paired ttest')
hold on
if p_post1_post2_LM1<.05
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
if p_post2_post3_LM1<.05
    plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k')
end
if p_post1_post3_LM1<.05
    plot(xt([1 3]), [1 1]*max(yt)*1.3, '-k',  mean(xt([1 3])), max(yt)*1.35, '*k')
end
hold off

%%%%%%%%%%% Unlesioned hemisphere %%%%%%%%%%%
figure
violin([TBL.Q_base_post1_diff_mFR(IX_UnLesion_LM1) TBL.Q_base_post2_diff_mFR(IX_UnLesion_LM1) TBL.Q_base_post3_diff_mFR(IX_UnLesion_LM1)], ...
    'xlabel',{'Post LDOPA','Peak 80','Post 60 min'})

