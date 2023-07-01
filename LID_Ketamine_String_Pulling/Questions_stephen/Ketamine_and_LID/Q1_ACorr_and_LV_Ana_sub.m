%% Corr and autocorr analysis.

ACb = AC_base_pre_early_late(:,:,1);
ACp1 = AC_base_pre_early_late(:,:,2);
ACp2 = AC_base_pre_early_late(:,:,3);
ACp3 = AC_base_pre_early_late(:,:,4);

 
% CORRELATIONS: Hypothesis: ketamine shifts cell-pair correlations and
% makes them weaker.
lbls = {'Base' 'Post 1' 'Post 2'};
figure
IXminp = min(ALL_R_p,[],2) < 0.01;
% IX = ALL_R_p(:,1) < 0.01;
subplot(2,1,1)
histogram_cowen({ALL_R(IXminp,1) ALL_R(IXminp,2) ALL_R(IXminp,3)},.005)
legend(lbls)
subplot(2,1,2)
histogram_cowen({ALL_R(:,2) - ALL_R(:,1) ALL_R(:,3) - ALL_R(:,1)  },.005,[],[],'pdf')
plot_vert_line_at_zero

% Hypothsis: Ketmaine makes positive correlations more positive and negative correlations more negative.
% it puts the net into a stronger state.
% I think this analysis suffers from regression to the mean if you are not
% careful.
%
% look at positive correlatins during pre
% sig correlations that were positive during pre. Did they go up or down in
% post.
IX = ALL_R_p(:,1) < 0.01 ;
% IXneg = ALL_R_p(:,1) < 0.01 & ALL_R(:,1)< 0;
IXpos = ALL_R(:,1)> 0;
IXneg = ALL_R(:,1)< 0;
figure
subplot(2,1,1)
histogram_cowen({ALL_R(:,2) - ALL_R(:,1) ALL_R(IX,2) - ALL_R(IX,1)  },.005,[],[],'pdf')
plot_vert_line_at_zero
subplot(2,1,2)
histogram_cowen({ALL_R(IXpos,2) - ALL_R(IXpos,1) ALL_R(IXneg,2) - ALL_R(IXneg,1) },.005,[],[],'pdf')
plot_vert_line_at_zero
[~,ptp] = ttest(ALL_R(IXpos,2) - ALL_R(IXpos,1));
[~,ptn] = ttest(ALL_R(IXneg,2) - ALL_R(IXneg,1));

% 
IXminp = min(ALL_R_p(:,1:2),[],2) < 0.01;

figure
subplot(2,1,1)
histogram(abs(ALL_R(IXminp,2)) - abs(ALL_R(IXminp,1)),60)
plot_vert_line_at_zero
pubify_figure_axis


xlabel('r')
subplot(2,1,2)
cnts = sum(ALL_R_p<0.01);
title('Difference r(sig) between post-pre')
bar( sum(ALL_R_p<0.01),'k')
pch = chi2_test(cnts, [mean(cnts) mean(cnts) mean(cnts) mean(cnts)])
ylabel('# significant correlations')
set(gca,'XTickLabel',lbls)
pubify_figure_axis

%% Autocorr.
mn_b = sum(AllQ{1},1)/(diff(Dset.intervals_around_evt_min(1,:))*60);
mn_p1 = sum(AllQ{2},1)/(diff(Dset.intervals_around_evt_min(2,:))*60);
mn_p2 = sum(AllQ{3},1)/(diff(Dset.intervals_around_evt_min(3,:))*60);

Error_bars([mn_b(:) mn_p1(:) mn_p2(:)])
kruskalwallis([mn_b(:) mn_p1(:) mn_p2(:)])

BASEIX = Dset.AC_x_ms > 120;
XIX = Dset.AC_x_ms < 60;
% Show each one.
if 0
    
    % Plot one at a time and PAUSE - press space to continue 
for ii = 1:Rows(ACb)
        figure
        plot(Dset.AC_x_ms(XIX),ACb(ii,XIX),'Color',cm(1,:),'LineWidth',4)
        hold on
        plot(Dset.AC_x_ms(XIX),ACp1(ii,XIX),'Color',cm(2,:),'LineWidth',4)
        plot(Dset.AC_x_ms(XIX),ACp2(ii,XIX),'Color',cm(3,:),'LineWidth',4)
        pubify_figure_axis
        legend(LBLS)
        legend boxoff
        pause
        close
        xlabel('ms')  
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot all autocorrs for all cells before and after ketamine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only look at autocorrs from cells that actually fire somewhat
GIXACR = min([mn_b(GIXright); mn_p1(GIXright); mn_p2(GIXright)])>0.1;
new_ACb = ACb(GIXACR,XIX);
new_ACp1 = ACp1(GIXACR,XIX);
new_ACp2 = ACp2(GIXACR,XIX);
% norm
base = mean(ACb(GIXACR,BASEIX),2);
% base = mean(ACb(GIX,XIX),2);
new_ACb = new_ACb - base;
new_ACp1 =new_ACp1 - base;
new_ACp2 =new_ACp2 - base;

new_ACb = new_ACb - mean(new_ACb,2);
new_ACp1 =new_ACp1 - mean(new_ACp1,2);
new_ACp2 =new_ACp2 - mean(new_ACp2,2);


[~,six] = sort(nanmax(new_ACb,[],2));

figure
subplot(1,3,1)
imagesc(Dset.AC_x_ms(XIX),[],new_ACb(six,:))
xlabel('ms')
ylabel('Neuron (sorted by peak)')
title('Baseline Norm Autocorr')
pubify_figure_axis
subplot(1,3,2)
imagesc(Dset.AC_x_ms(XIX),[],new_ACp1(six,:))
pubify_figure_axis
title('Post 1')
subplot(1,3,3)
imagesc(Dset.AC_x_ms(XIX),[],new_ACp2(six,:))
pubify_figure_axis
title('Post 2')

equalize_color_axes
colorbar_label
sgtitle('Right hemisphere LDo&ket SHAM')

GIXACL = min([mn_b(GIXleft); mn_p1(GIXleft); mn_p2(GIXleft)])>0.1;
new_ACb = ACb(GIXACL,XIX);
new_ACp1 = ACp1(GIXACL,XIX);
new_ACp2 = ACp2(GIXACL,XIX);
% norm
base = mean(ACb(GIXACL,BASEIX),2);
% base = mean(ACb(GIX,XIX),2);
new_ACb = new_ACb - base;
new_ACp1 =new_ACp1 - base;
new_ACp2 =new_ACp2 - base;

new_ACb = new_ACb - mean(new_ACb,2);
new_ACp1 =new_ACp1 - mean(new_ACp1,2);
new_ACp2 =new_ACp2 - mean(new_ACp2,2);


[~,six] = sort(nanmax(new_ACb,[],2));

figure
subplot(1,3,1)
imagesc(Dset.AC_x_ms(XIX),[],new_ACb(six,:))
xlabel('ms')
ylabel('Neuron (sorted by peak)')
title('Baseline Norm Autocorr')
pubify_figure_axis
subplot(1,3,2)
imagesc(Dset.AC_x_ms(XIX),[],new_ACp1(six,:))
pubify_figure_axis
title('Post 1')
subplot(1,3,3)
imagesc(Dset.AC_x_ms(XIX),[],new_ACp2(six,:))
pubify_figure_axis
title('Post 2')

equalize_color_axes
colorbar_label
sgtitle('left hemisphere LDo&ket SHAM')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot_confidence_intervals(Dset.AC_x_ms(XIX),new_ACb,[],cm(1,:))
hold on
plot_confidence_intervals(Dset.AC_x_ms(XIX),new_ACp1,[],cm(2,:))
plot_confidence_intervals(Dset.AC_x_ms(XIX),new_ACp2,[],cm(3,:))
legend_color_text(lbls,cm)
pubify_figure_axis
xlabel('ms')
ylabel('Normalized Rate')
title('Mean Auto Corrs Right hemisphere')

figure
plot_confidence_intervals(Dset.AC_x_ms(XIX),new_ACb,[],cm(1,:))
hold on
plot_confidence_intervals(Dset.AC_x_ms(XIX),new_ACp1,[],cm(2,:))
plot_confidence_intervals(Dset.AC_x_ms(XIX),new_ACp2,[],cm(3,:))
legend_color_text(lbls,cm)
pubify_figure_axis
xlabel('ms')
ylabel('Normalized Rate')
title('Mean Auto Corrs Left hemisphere')
% relaated is local variance
clrs = lines(5);

figure
subplot(2,1,1)
p1 = signrank(TBL.LocVar_post1(GIXACR) -TBL.LocVar_base(GIXACR));
mn1 = nanmean(TBL.LocVar_post1(GIXACR) -TBL.LocVar_base(GIXACR));
p2 = signrank(TBL.LocVar_post2(GIXACR) -TBL.LocVar_base(GIXACR));
mn2 = nanmean(TBL.LocVar_post2(GIXACR) -TBL.LocVar_base(GIXACR));
p2p1 = signrank(TBL.LocVar_post2(GIXACR) -TBL.LocVar_post1(GIXACR));
mn2m1 = nanmean(TBL.LocVar_post2(GIXACR) -TBL.LocVar_post1(GIXACR));

histogram_cowen({TBL.LocVar_base(GIXACR) TBL.LocVar_post1(GIXACR) TBL.LocVar_post2(GIXACR)},.1 )
title(sprintf('Right hemisphere m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f,  m2m1=%1.2f, p2p1=%1.3f',mn1,p1,mn2,p2,mn2m1,p2mp1))
plot_ref_line(nanmean(TBL.LocVar_base(GIXACR)),'line_width',2,'style','-','color',clrs(1,:))
plot_ref_line(nanmean(TBL.LocVar_post1(GIXACR)),'line_width',2,'style','-','color',clrs(2,:))
plot_ref_line(nanmean(TBL.LocVar_post2(GIXACR)),'line_width',2,'style','-','color',clrs(3,:))
legend('Baseline','L-DOPA','Ketamine');legend boxoff
pubify_figure_axis

subplot(2,1,2)
p1 = signrank(TBL.LocVar_post1(GIXACL) -TBL.LocVar_base(GIXACL));
mn1 = nanmean(TBL.LocVar_post1(GIXACL) -TBL.LocVar_base(GIXACL));
p2 = signrank(TBL.LocVar_post2(GIXACL) -TBL.LocVar_base(GIXACL));
mn2 = nanmean(TBL.LocVar_post2(GIXACL) -TBL.LocVar_base(GIXACL));
p2p1 = signrank(TBL.LocVar_post2(GIXACL) -TBL.LocVar_post1(GIXACL));
mn2m1 = nanmean(TBL.LocVar_post2(GIXACL) -TBL.LocVar_post1(GIXACL));

histogram_cowen({TBL.LocVar_base(GIXACL) TBL.LocVar_post1(GIXACL) TBL.LocVar_post2(GIXACL)},.1 )
title(sprintf('Left hemisphere m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f, m2m1=%1.2f, p2p1=%1.3f',mn1,p1,mn2,p2,mn2m1,p2mp1))
xlabel('Local Variance')
plot_ref_line(nanmean(TBL.LocVar_base(GIXACL)),'line_width',2,'style','-','color',clrs(1,:))
plot_ref_line(nanmean(TBL.LocVar_post1(GIXACL)),'line_width',2,'style','-','color',clrs(2,:))
plot_ref_line(nanmean(TBL.LocVar_post2(GIXACL)),'line_width',2,'style','-','color',clrs(3,:))
legend('Baseline','L-DOPA','Ketamine');legend boxoff
pubify_figure_axis


figure
subplot(2,1,1)
histogram_cowen({(TBL.LocVar_post1(GIXACR) -TBL.LocVar_base(GIXACR)) (TBL.LocVar_post2(GIXACR) -TBL.LocVar_base(GIXACR)) ...
    (TBL.LocVar_post2(GIXACR) -TBL.LocVar_post1(GIXACR)) },.05 )
plot_vert_line_at_zero
legend('L-DOPA-Baseline','Ketamine-Baseline', 'L-DOPA-Ketamine');legend boxoff
% plot_vert_line_at_zero
% plot_vert_line_at_zero(nanmean(TBL.LocVar_post2(GIXright) -TBL.LocVar_base(GIXright)))
% plot_vert_line_at_zero(nanmean(TBL.LocVar_post1(GIXright) -TBL.LocVar_base(GIXright)))
% plot_vert_line_at_zero(nanmean(TBL.LocVar_post2(GIXright) -TBL.LocVar_post1(GIXright)))
xlabel('Difference measure of Local Variance')
title('Lesioned hemisphere')
pubify_figure_axis
signrank((TBL.LocVar_post2(GIXACR) -TBL.LocVar_base(GIXACR)))
signrank((TBL.LocVar_post1(GIXACR) -TBL.LocVar_base(GIXACR)))
signrank((TBL.LocVar_post2(GIXACR) -TBL.LocVar_post1(GIXACR)))
signrank((TBL.LocVar_post3(GIXACR) -TBL.LocVar_post1(GIXACR)))

subplot(2,1,2)
histogram_cowen({(TBL.LocVar_post1(GIXACL) -TBL.LocVar_base(GIXACL)) (TBL.LocVar_post2(GIXACL) -TBL.LocVar_base(GIXACL)) ...
    (TBL.LocVar_post2(GIXACL) -TBL.LocVar_post1(GIXACL))},.05 )
plot_vert_line_at_zero
legend('L-DOPA-Baseline','Ketamine-Baseline', 'L-DOPA-Ketamine');legend boxoff
xlabel('Difference measure of Local Variance')
title('Un-Lesioned hemisphere')
pubify_figure_axis
signrank((TBL.LocVar_post2(GIXACL) -TBL.LocVar_base(GIXACL)))
signrank((TBL.LocVar_post1(GIXACL) -TBL.LocVar_base(GIXACL)))
signrank((TBL.LocVar_post2(GIXACL) -TBL.LocVar_post1(GIXACL)))
signrank((TBL.LocVar_post3(GIXACL) -TBL.LocVar_post1(GIXACL)))