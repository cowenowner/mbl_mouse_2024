
%Control Motor Cortex
IX_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control';
%Control Striatum
IX_CS = categorical(TBL.BrainRegion) == 'Striatum' & categorical(TBL.Group) == 'Control' & TBL.Frate_base > 0 & TBL.Frate_post1 > 0;
%LID
IX_Lesion_LM1 =  categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Hemisphere == 'R' ; % 55 neurons 
IX_UnLesion_LM1 =  categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Hemisphere == 'L'; % 55 neurons 

%LID Striatum
IX_LS =  categorical(TBL.BrainRegion) == 'Striatum' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Frate_base > 0 & TBL.Frate_post1 > 0; % 55 neurons 

TBL.m_LocVar = nanmean([TBL.LocVar_base TBL.LocVar_post2],2);
TBL.m_Frate = nanmean([TBL.Frate_base TBL.Frate_post2],2);
figure
subplot(411)
scatter(TBL.m_Frate(IX_CM1), TBL.Frate_post1mbase(IX_CM1))
lsline
xlabel('mean FR')
ylabel('FR post1 - base')
title('FR change from baseline in Control M1')
plot_ref_line(0,'orientation', 'horiz')

subplot(412)
scatter(TBL.m_LocVar(IX_CM1), TBL.LocVar_post1mbase(IX_CM1))
xlabel('mean Local Variance')
ylabel('Loc Var post1 - base')
title('Local Variance change from Loc Var baseline in Control M1')
lsline
plot_ref_line(0,'orientation', 'horiz')

subplot(413)
scatter(TBL.m_LocVar(IX_CM1), TBL.Frate_post1mbase(IX_CM1))
xlabel('mean Local Variance')
ylabel('FR post1 - base')
title('FR change from Loc Var baseline in Control M1')
lsline
plot_ref_line(0,'orientation', 'horiz')

subplot(414)
scatter( TBL.LocVar_post1mbase(IX_CM1),  TBL.Frate_post1mbase(IX_CM1))
xlabel('Loc Var post1 - base')
ylabel('FR post1 - base')
title('FR change from Loc Var change in Control M1')
lsline
plot_ref_line(0,'orientation', 'horiz')

%% Base and post plotted for FR and Loc Var
% Control M1
figure
subplot(121)
scatter(TBL.Frate_base(IX_CM1 & IR_resp), TBL.LocVar_base(IX_CM1 & IR_resp))
lsline
ylabel('Loc Var base')
xlabel('FR base')
title('baseline in Control M1')

subplot(122)
scatter(TBL.Frate_post1(IX_CM1 & IR_resp), TBL.LocVar_post1(IX_CM1 & IR_resp))
lsline
ylabel('Loc Var post1')
xlabel('FR post1')
title('Post ketamine in Control M1')

% M1 %
% plot the histograms for the firing rat and loc var
figure
% subplot(211)
histogram_cowen({TBL.LocVar_salmbase(IX_CM1 & IR_resp) TBL.LocVar_postket1mbase(IX_CM1 & IR_resp) ...
    TBL.LocVar_postket2mbase(IX_CM1 & IR_resp)},.05)
xlabel('Loc Var')
title('LocVar CM1')
pubify_figure_axis
ax=gca;
ax.FontSize=35;
% set(gca,'xlim',[-60 120])
plot_ref_line(0,'line_width',2,'style',':')

subplot(212)
histogram_cowen({TBL.Frate_salmbase(IX_CM1 & IR_resp) TBL.Frate_postket1mbase(IX_CM1 & IR_resp) ...
    TBL.Frate_postket2mbase(IX_CM1 & IR_resp)},1)
xlabel('FR post sal and post ket1 and post ket2')
title('Firing rate change post sal and psot ket relative to baseline Control M1')
plot_ref_line(0)

% Test if the Local Variance is different
[h,p_CM1(1)] = ttest(TBL.LocVar_salmbase(IX_CM1 & IR_resp));
[h,p_CM1(2)] = ttest(TBL.LocVar_postket1mbase(IX_CM1 & IR_resp));
[h,p_CM1(3)] = ttest(TBL.LocVar_postket2mbase(IX_CM1 & IR_resp));

d(1) = Cohens_d(TBL.LocVar_salmbase(IX_CM1 & IR_resp));
d(2) = Cohens_d(TBL.LocVar_postket1mbase(IX_CM1 & IR_resp));
d(3) = Cohens_d(TBL.LocVar_postket2mbase(IX_CM1 & IR_resp));

p_CM1_comp(1) = ranksum(TBL.LocVar_postket1mbase(IX_CM1 & IR_resp),TBL.LocVar_salmbase(IX_CM1 & IR_resp));
p_CM1_comp(2) = ranksum(TBL.LocVar_postket1mbase(IX_CM1 & IR_resp),TBL.LocVar_postket2mbase(IX_CM1 & IR_resp));
p_CM1_comp(3) = ranksum(TBL.LocVar_postket2mbase(IX_CM1 & IR_resp),TBL.LocVar_salmbase(IX_CM1 & IR_resp));

d_comp(1) = Cohens_d(TBL.LocVar_postket1mbase(IX_CM1 & IR_resp),TBL.LocVar_salmbase(IX_CM1 & IR_resp));
d_comp(2) = Cohens_d(TBL.LocVar_postket2mbase(IX_CM1 & IR_resp),TBL.LocVar_postket1mbase(IX_CM1 & IR_resp));
d_comp(3) = Cohens_d(TBL.LocVar_postket2mbase(IX_CM1 & IR_resp),TBL.LocVar_salmbase(IX_CM1 & IR_resp));


%%% Control Straitum %%
figure
subplot(121)
scatter(TBL.Frate_base(IX_CS & IR_resp), TBL.LocVar_base(IX_CS & IR_resp))
lsline
ylabel('Loc Var base')
xlabel('FR base')
title('baseline in Control Striatum')

subplot(122)
scatter(TBL.Frate_post1(IX_CS & IR_resp), TBL.LocVar_post1(IX_CS & IR_resp))
lsline
ylabel('Loc Var post1')
xlabel('FR post1')
title('Post ketamine in Control Striatum')

% M1 %
% plot the histograms for the firing rat and loc var
figure
% subplot(211)
histogram_cowen({TBL.LocVar_salmbase(IX_CS & IR_resp) TBL.LocVar_postket1mbase(IX_CS & IR_resp) ...
    TBL.LocVar_postket2mbase(IX_CS & IR_resp)},.05)
xlabel('Loc Var')
title('LocVar Control Striatum')
pubify_figure_axis
ax=gca;
ax.FontSize=35;
% set(gca,'xlim',[-60 120])
plot_ref_line(0,'line_width',2,'style',':')

subplot(212)
histogram_cowen({TBL.Frate_salmbase(IX_CS & IR_resp) TBL.Frate_postket1mbase(IX_CS & IR_resp) ...
    TBL.Frate_postket2mbase(IX_CS & IR_resp)},1)
xlabel('FR post sal and post ket1 and post ket2')
title('Firing rate change post sal and psot ket relative to baseline Control Striatum')
plot_ref_line(0)

% Test if the Local Variance is different
[h,p_CS(1)] = ttest(TBL.LocVar_salmbase(IX_CS & IR_resp));
[h,p_CS(2)] = ttest(TBL.LocVar_postket1mbase(IX_CS & IR_resp));
[h,p_CS(3)] = ttest(TBL.LocVar_postket2mbase(IX_CS & IR_resp));

d(1) = Cohens_d(TBL.LocVar_salmbase(IX_CS & IR_resp));
d(2) = Cohens_d(TBL.LocVar_postket1mbase(IX_CS & IR_resp));
d(3) = Cohens_d(TBL.LocVar_postket2mbase(IX_CS & IR_resp));

p_CS_comp(1) = ranksum(TBL.LocVar_postket1mbase(IX_CS & IR_resp),TBL.LocVar_salmbase(IX_CS & IR_resp));
p_CS_comp(2) = ranksum(TBL.LocVar_postket1mbase(IX_CS & IR_resp),TBL.LocVar_postket2mbase(IX_CS & IR_resp));
p_CS_comp(3) = ranksum(TBL.LocVar_postket2mbase(IX_CS & IR_resp),TBL.LocVar_salmbase(IX_CS & IR_resp));

d_comp(1) = Cohens_d(TBL.LocVar_postket1mbase(IX_CS & IR_resp),TBL.LocVar_salmbase(IX_CS & IR_resp));
d_comp(2) = Cohens_d(TBL.LocVar_postket2mbase(IX_CS & IR_resp),TBL.LocVar_postket1mbase(IX_CS & IR_resp));
d_comp(3) = Cohens_d(TBL.LocVar_postket2mbase(IX_CS & IR_resp),TBL.LocVar_salmbase(IX_CS & IR_resp));

% Lesioned hemisphere tests
[h,p_Lesion_LM1(1)] = ttest(TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));
[h,p_Lesion_LM1(2)] = ttest(TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp));
[h,p_Lesion_LM1(3)] = ttest(TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp));

d(1) = Cohens_d(TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));
d(2) = Cohens_d(TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp));
d(3) = Cohens_d(TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp));

p_Lesion_LM1_comp(1) = ranksum(TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));
p_Lesion_LM1_comp(2) = ranksum(TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp));
p_Lesion_LM1_comp(3) = ranksum(TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));

d_comp(1) = Cohens_d(TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));
d_comp(2) = Cohens_d(TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp));
d_comp(3) = Cohens_d(TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));

% Unlesioned hemisphere tests
[h,p_UnLesion_LM1(1)] = ttest(TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLesion_LM1(2)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLesion_LM1(3)] = ttest(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));

d(1) = Cohens_d(TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
d(2) = Cohens_d(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp));
d(3) = Cohens_d(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));

p_UnLesion_LM1_comp(1) = ranksum(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
p_UnLesion_LM1_comp(2) = ranksum(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));
p_UnLesion_LM1_comp(3) = ranksum(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));

d_comp(1) = Cohens_d(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
d_comp(2) = Cohens_d(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp));
d_comp(3) = Cohens_d(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));

% 
% figure
% histogram_cowen({TBL.Frate_base TBL.Frate_post1},.05)
% xlabel('FR pre and post')
% title('Firing rate M1 Control')

% cohen's d comparison of Loc Var and FR change
LocVar_cohens_d = Cohens_d(Locvar_post1base_CM1);
FR_cohens_d = Cohens_d(FR_post1base_CM1);
% cohen's d comparison of Loc Var and FR change
LocVar_cohens_d_CS = Cohens_d(Locvar_post1base_CS);
FR_cohens_d_CS = Cohens_d(FR_post1base_CS);
% Correlation
[rho,pval] = corr([TBL.Frate_post1(IX_CS) TBL.LocVar_post1(IX_CS)],'rows','complete');

% rank sum test for the post-pre histogram 
% loc var control post ket
IX_ket_varBzero_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.LocVar_postket1mbase < 0;
IX_ket_varAzero_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.LocVar_postket1mbase > 0;

p_var_ket_CM1 = ranksum(TBL.LocVar_postket1mbase(IX_ket_varBzero_CM1),TBL.LocVar_postket1mbase(IX_ket_varAzero_CM1));
% loc var post saline
IX_sal_varBzero_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.LocVar_salmbase < 0;
IX_sal_varAzero_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.LocVar_salmbase > 0;

p_var_sal_CM1 = ranksum(TBL.LocVar_salmbase(IX_sal_varBzero_CM1),TBL.LocVar_salmbase(IX_sal_varAzero_CM1));
% loc var LID
IX_varBzero_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Frate_base > 0 ...
    & TBL.Frate_post2 > 0 & TBL.LocVar_post2mbase < 0;
IX_varAzero_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Frate_base > 0 ...
    & TBL.Frate_post2 > 0 & TBL.LocVar_post2mbase > 0;

p_var_LM1 = ranksum(TBL.LocVar_post1mbase(IX_varBzero_LM1),TBL.LocVar_post1mbase(IX_varAzero_LM1));
% f rate change control
IX_frateBzero_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Frate_base > 0 ...
    & TBL.Frate_post1 > 0 & TBL.Frate_post1mbase < 0;
IX_frateAzero_CM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == 'Control' & TBL.Frate_base > 0 ...
    & TBL.Frate_post1 > 0 & TBL.Frate_post1mbase > 0;

p_frate_CM1 = ranksum(TBL.Frate_post1mbase(IX_frateBzero_CM1),TBL.Frate_post1mbase(IX_frateAzero_CM1));
%frate change LID
IX_frateBzero_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Frate_base > 0 ...
    & TBL.Frate_post2 > 0 & TBL.Frate_post2mbase < 0;
IX_frateAzero_LM1 = categorical(TBL.BrainRegion) == 'M1' & categorical(TBL.Group) == '6ODHA_LID' & TBL.Frate_base > 0 ...
    & TBL.Frate_post2 > 0 & TBL.Frate_post2mbase > 0;

p_frate_LM1 = ranksum(TBL.Frate_post1mbase(IX_frateBzero_LM1),TBL.Frate_post1mbase(IX_frateAzero_LM1));

% % Striatum %
% % plot the histograms for the firing rate and loc var
% figure
% subplot(211)
% histogram_cowen({Locvar_post1base_CS},.05)
% xlabel('Loc Var post - base')
% title('LocVar post1mbase Striatum Control')
% 
% subplot(212)
% histogram_cowen({FR_post1base_CS},.05)
% xlabel('FR post - base')
% title('Firing rate change Striatum Control')


% M1 and Striatum %
% hist ISI
figure
subplot(221)
HistISI(Dset.SP(3).t_uS/100)
title('ISI NRN 3 Session 12 Control M1')
t_uS_diff = diff(Dset.SP(24).t_uS);
subplot(222)
plot(t_uS_diff)
title('Diff NRN 3 Session 12 Control M1')
subplot(223)
plot(t_uS_diff)
title('Diff NRN 8 Session 12 Control M1')

figure 
subplot(221)
plot(t_uS_diff)
title('Diff NRN 24 Session 12 Control M1')
subplot(222)
histogram(t_uS_diff/1000,10)
title('Diff NRN 24 Session 12 Control M1')
% Auto corr
[AC,x] = AutoCorr(Dset.SP(3).t_uS/1000,5,300);
subplot(222)
plot(x,AC)
title('AutoCorr NRN 3 Session 12 Control M1')
% hist ISI
figure
subplot (221)
HistISI(Dset.SP(24).t_uS/100)
title('ISI NRN 24 Session 12 Control m1')
% Auto corr
[AC,x] = AutoCorr(Dset.SP(19).t_uS/1000,5,300);
subplot(222)
plot(x,AC)
title('AutoCorr NRN 19 Session 12 Control Striatum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LID animal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M1

figure
subplot(411)
scatter(TBL.m_Frate(IX_LM1), TBL.Frate_post1mbase(IX_LM1))
lsline
xlabel('mean FR')
ylabel('FR post1 - base')
title('FR change from mean in LID M1')
plot_ref_line(0,'orientation', 'horiz')

subplot(412)
scatter(TBL.m_LocVar(IX_LM1), TBL.LocVar_post1mbase(IX_LM1))
xlabel('mean Local Variance')
ylabel('Loc Var post1 - base')
title('Local Variance change from mean in LID M1')
lsline
plot_ref_line(0,'orientation', 'horiz')

subplot(413)
scatter(TBL.m_LocVar(IX_LM1), TBL.Frate_post1mbase(IX_LM1))
xlabel('mean Local Variance')
ylabel('FR post1 - base')
title('FR change from Loc Var mean in LID M1')
lsline
plot_ref_line(0,'orientation', 'horiz')

subplot(414)
scatter(TBL.LocVar_post1mbase(IX_LM1), TBL.Frate_post1mbase(IX_LM1))
xlabel('Local Variance post - base')
ylabel('FR post1 - base')
title('FR change from Loc Var change in LID M1')
lsline
plot_ref_line(0,'orientation', 'horiz')

% Hist ISI
figure
subplot(221)
HistISI(Dset.SP(20).t_uS/100)
title('ISI NRN 20 Session 9 LID M1 Left Hemisphere')
% Auto corr
[AC,x] = AutoCorr(Dset.SP(20).t_uS/1000,5,300);
subplot(222)
plot(x,AC)
title('AutoCorr NRN 20 Session 9 LID M1 Left Hemisphere')
% Hist ISI
subplot(223)
HistISI(Dset.SP(17).t_uS/100)
title('ISI NRN 17 Session 16 LID M1 Right Hemisphere')
% Auto corr
[AC,x] = AutoCorr(Dset.SP(17).t_uS/1000,5,300);
subplot(224)
plot(x,AC)
title('AutoCorr NRN 17 Session 16 LID M1 Right Hemisphere')

% plot the histograms for the firing rat and loc var
% Pre and Post KETAMINE %%%%
%%%%%%lesioned hem %%%%%%%%
figure
% subplot(211)
histogram_cowen({TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp),...
    TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp)},.05)
xlabel('Loc Var')
title('LocVar LM1 Lesion')
pubify_figure_axis
ax=gca;
ax.FontSize=35;
% set(gca,'xlim',[-60 120])
plot_ref_line(0,'line_width',2,'style',':')

subplot(212)
histogram_cowen({TBL.Frate_salmbase(IX_Lesion_LM1 & IR_resp),TBL.Frate_postket1mbase(IX_Lesion_LM1 & IR_resp),...
    TBL.Frate_postket2mbase(IX_Lesion_LM1 & IR_resp)},1)
xlabel('FR postsal psotket and postket2 relative to baseline')
title('Firing rate change post saline and ketamine M1 LID Lesion')
plot_ref_line(0)

% ttest 
[h,p_Les_LM1(1)] = ttest(TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));
[h,p_Les_LM1(2)] = ttest(TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp));
[h,p_Les_LM1(3)] = ttest(TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp));

[h,p_Les_LM1_comp(1),ci,stats(1)] = ttest(TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));
[h,p_Les_LM1_comp(2),ci,stats(2)] = ttest(TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp));
[h,p_Les_LM1_comp(3),ci,stats(3)] = ttest(TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp));

%%%%% Unlesioned hem %%%%%%
figure
% subplot(211)
histogram_cowen({TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),...
    TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp)},.05)
xlabel('Loc Var')
title('LocVar LM1 UnLesion')
pubify_figure_axis
ax=gca;
ax.FontSize=25;
% set(gca,'xlim',[-60 120])
plot_ref_line(0,'line_width',2,'style',':')

subplot(212)
histogram_cowen({TBL.Frate_salmbase(IX_UnLesion_LM1 & IR_resp),TBL.Frate_postket1mbase(IX_UnLesion_LM1 & IR_resp),...
    TBL.Frate_postket2mbase(IX_UnLesion_LM1 & IR_resp)},1)
xlabel('FR postsal psotket and postket2 relative to baseline')
title('Firing rate change post saline and ketamine M1 LID UnLesion')
plot_ref_line(0)

% ttest
[h,p_UnLes_LM1(1)] = ttest(TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLes_LM1(2)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLes_LM1(3)] = ttest(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));

[h,p_UnLes_LM1_comp(1),ci,stats(1)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLes_LM1_comp(2),ci,stats(2)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLes_LM1_comp(3),ci,stats(3)] = ttest(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));


% baseline and post of FR and LocVar
%%% lesioned %%%%
figure
subplot(121)
scatter(TBL.Frate_base(IX_Lesion_LM1 & IR_resp), TBL.LocVar_base(IX_Lesion_LM1 & IR_resp))
lsline
ylabel('Loc Var base')
xlabel('FR base')
title('baseline in lesioned LID M1')

subplot(122)
scatter(TBL.Frate_post1(IX_Lesion_LM1 & IR_resp), TBL.LocVar_post1(IX_Lesion_LM1 & IR_resp))
lsline
ylabel('Loc Var post1')
xlabel('FR post1')
title('Post ket in lesioned LID M1')

[rho,pval] = corr([TBL.Frate_base(IX_LM1) TBL.LocVar_base(IX_LM1)],'rows','complete');

%%% Unlesioned %%%%%%
figure
subplot(121)
scatter(TBL.Frate_base(IX_UnLesion_LM1), TBL.LocVar_base(IX_UnLesion_LM1))
lsline
ylabel('Loc Var base')
xlabel('FR base')
title('baseline in Unlesioned LID M1')

subplot(122)
scatter(TBL.Frate_post1(IX_UnLesion_LM1), TBL.LocVar_post1(IX_UnLesion_LM1))
lsline
ylabel('Loc Var post1')
xlabel('FR post1')
title('Post ket in Unlesioned LID M1')

%% % Pre and Post LDOPA %%%%
% plot the histograms for the firing rat and loc var
%%%%%%lesioned hem %%%%%%%%
figure
% subplot(211)
histogram_cowen({TBL.LocVar_salmbase(IX_Lesion_LM1 & IR_resp),TBL.LocVar_postket1mbase(IX_Lesion_LM1 & IR_resp),...
    TBL.LocVar_postket2mbase(IX_Lesion_LM1 & IR_resp)},.05)
xlabel('Loc Var')
title('LocVar LM1 Lesion')
pubify_figure_axis
ax=gca;
ax.FontSize=25;
% set(gca,'xlim',[-60 120])
plot_ref_line(0,'line_width',2,'style',':')

subplot(212)
histogram_cowen({TBL.Frate_salmbase(IX_Lesion_LM1 & IR_resp),TBL.Frate_postket1mbase(IX_Lesion_LM1 & IR_resp),...
    TBL.Frate_postket2mbase(IX_Lesion_LM1 & IR_resp)},1)
xlabel('FR postLDOPA peak80 postket relative to baseline')
title('Firing rate change post LDOPA, peak80 and ketamine M1 LID Lesion')
plot_ref_line(0)
%ttest
[h,p_Les_LDO_LM1(1)] = ttest(TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
[h,p_Les_LDO_LM1(2)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp));
[h,p_Les_LDO_LM1(3)] = ttest(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));

[h,p_Les_LDO_LM1_comp(1),ci,stats(1)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
[h,p_Les_LDO_LM1_comp(2),ci,stats(2)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));
[h,p_Les_LDO_LM1_comp(3),ci,stats(3)] = ttest(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));


%%%%% Unlesioned hem %%%%%%
figure
% subplot(211)
histogram_cowen({TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),...
    TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp)},.05)
xlabel('Loc Var')
title('LocVar LM1 UnLesion')
pubify_figure_axis
ax=gca;
ax.FontSize=25;
% set(gca,'xlim',[-60 120])
plot_ref_line(0,'line_width',2,'style',':')

subplot(212)
histogram_cowen({TBL.Frate_salmbase(IX_UnLesion_LM1 & IR_resp),TBL.Frate_postket1mbase(IX_UnLesion_LM1 & IR_resp),...
    TBL.Frate_postket2mbase(IX_UnLesion_LM1 & IR_resp)},1)
xlabel('FR postLDOPA peak80 postket relative to baseline')
title('Firing rate change post LDOPA, peak80 and ketamine M1 LID UnLesion')
plot_ref_line(0)
%ttest
[h,p_UnLes_LDO_LM1(1)] = ttest(TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLes_LDO_LM1(2)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLes_LDO_LM1(3)] = ttest(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));

[h,p_UnLes_LDO_LM1_comp(1),ci,stats(1)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLes_LDO_LM1_comp(2),ci,stats(2)] = ttest(TBL.LocVar_postket1mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp));
[h,p_UnLes_LDO_LM1_comp(3),ci,stats(3)] = ttest(TBL.LocVar_postket2mbase(IX_UnLesion_LM1 & IR_resp),TBL.LocVar_salmbase(IX_UnLesion_LM1 & IR_resp));

% baseline and post of FR and LocVar
%%%%% lesioned hemisphere %%%%%%%
figure
subplot(131)
scatter(TBL.Frate_base(IX_Lesion_LM1), TBL.LocVar_base(IX_Lesion_LM1))
lsline
ylabel('Loc Var base')
xlabel('FR base')
title('baseline in Lesioned LID M1')

subplot(132)
scatter(TBL.Frate_post2(IX_Lesion_LM1), TBL.LocVar_post2(IX_Lesion_LM1))
lsline
ylabel('Loc Var peak 80')
xlabel('FR peak80')
title('Peak 80 Lesioned in LID M1')

subplot(133)
scatter(TBL.Frate_post3(IX_Lesion_LM1), TBL.LocVar_post3(IX_Lesion_LM1))
lsline
ylabel('Loc Var post ket')
xlabel('FR post ket')
title('Post Ket Lesioned in LID M1')
%%% Unlesioend %%%%%%%
figure
subplot(131)
scatter(TBL.Frate_base(IX_UnLesion_LM1), TBL.LocVar_base(IX_UnLesion_LM1))
lsline
ylabel('Loc Var base')
xlabel('FR base')
title('baseline in UnLesioned LID M1')

subplot(132)
scatter(TBL.Frate_post2(IX_UnLesion_LM1), TBL.LocVar_post2(IX_UnLesion_LM1))
lsline
ylabel('Loc Var peak 80')
xlabel('FR peak80')
title('Peak 80 UnLesioned in LID M1')

subplot(133)
scatter(TBL.Frate_post3(IX_UnLesion_LM1), TBL.LocVar_post3(IX_UnLesion_LM1))
lsline
ylabel('Loc Var post ket')
xlabel('FR post ket')
title('Post Ket UnLesioned in LID M1')
