GIXrightPK = 1;
GIXleftPK = 1;

figure
subplot(2,1,1)
p1 = signrank(TBL.LocVar_post1(GIXrightPK) -TBL.LocVar_base(GIXrightPK));
mn1 = nanmean(TBL.LocVar_post1(GIXrightPK) -TBL.LocVar_base(GIXrightPK));
p2 = signrank(TBL.LocVar_post2(GIXrightPK) -TBL.LocVar_base(GIXrightPK));
mn2 = nanmean(TBL.LocVar_post2(GIXrightPK) -TBL.LocVar_base(GIXrightPK));
p2p1 = signrank(TBL.LocVar_post2(GIXrightPK) -TBL.LocVar_post1(GIXrightPK));
mn2m1 = nanmean(TBL.LocVar_post2(GIXrightPK) -TBL.LocVar_post1(GIXrightPK));

histogram_cowen({TBL.LocVar_base(GIXrightPK) TBL.LocVar_post1(GIXrightPK) TBL.LocVar_post2(GIXrightPK)},.1 )
title(sprintf('Right hemisphere m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f,  m2m1=%1.2f, p2p1=%1.3f',mn1,p1,mn2,p2,mn2m1,p2mp1))
plot_ref_line(nanmean(TBL.LocVar_base(GIXrightPK)),'line_width',2,'style','-','color',clrs(1,:))
plot_ref_line(nanmean(TBL.LocVar_post1(GIXrightPK)),'line_width',2,'style','-','color',clrs(2,:))
plot_ref_line(nanmean(TBL.LocVar_post2(GIXrightPK)),'line_width',2,'style','-','color',clrs(3,:))
legend('Baseline','L-DOPA','Ketamine');legend boxoff
pubify_figure_axis

subplot(2,1,2)
p1 = signrank(TBL.LocVar_post1(GIXleftPK) -TBL.LocVar_base(GIXleftPK));
mn1 = nanmean(TBL.LocVar_post1(GIXleftPK) -TBL.LocVar_base(GIXleftPK));
p2 = signrank(TBL.LocVar_post2(GIXleftPK) -TBL.LocVar_base(GIXleftPK));
mn2 = nanmean(TBL.LocVar_post2(GIXleftPK) -TBL.LocVar_base(GIXleftPK));
p2p1 = signrank(TBL.LocVar_post2(GIXleftPK) -TBL.LocVar_post1(GIXleftPK));
mn2m1 = nanmean(TBL.LocVar_post2(GIXleftPK) -TBL.LocVar_post1(GIXleftPK));

histogram_cowen({TBL.LocVar_base(GIXleftPK) TBL.LocVar_post1(GIXleftPK) TBL.LocVar_post2(GIXleftPK)},.1 )
title(sprintf('Left hemisphere m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f, m2m1=%1.2f, p2p1=%1.3f',mn1,p1,mn2,p2,mn2m1,p2mp1))
xlabel('Local Variance')
plot_ref_line(nanmean(TBL.LocVar_base(GIXleftPK)),'line_width',2,'style','-','color',clrs(1,:))
plot_ref_line(nanmean(TBL.LocVar_post1(GIXleftPK)),'line_width',2,'style','-','color',clrs(2,:))
plot_ref_line(nanmean(TBL.LocVar_post2(GIXleftPK)),'line_width',2,'style','-','color',clrs(3,:))
legend('Baseline','L-DOPA','Ketamine');legend boxoff
pubify_figure_axis

% IN

figure
subplot(2,1,1)
p1 = signrank(TBL.LocVar_post1(GIXrightIN) -TBL.LocVar_base(GIXrightIN));
mn1 = nanmean(TBL.LocVar_post1(GIXrightIN) -TBL.LocVar_base(GIXrightIN));
p2 = signrank(TBL.LocVar_post2(GIXrightIN) -TBL.LocVar_base(GIXrightIN));
mn2 = nanmean(TBL.LocVar_post2(GIXrightIN) -TBL.LocVar_base(GIXrightIN));
p2p1 = signrank(TBL.LocVar_post2(GIXrightIN) -TBL.LocVar_post1(GIXrightIN));
mn2m1 = nanmean(TBL.LocVar_post2(GIXrightIN) -TBL.LocVar_post1(GIXrightIN));

histogram_cowen({TBL.LocVar_base(GIXrightIN) TBL.LocVar_post1(GIXrightIN) TBL.LocVar_post2(GIXrightIN)},.1 )
title(sprintf('Right hemisphere m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f,  m2m1=%1.2f, p2p1=%1.3f',mn1,p1,mn2,p2,mn2m1,p2mp1))
plot_ref_line(nanmean(TBL.LocVar_base(GIXrightIN)),'line_width',2,'style','-','color',clrs(1,:))
plot_ref_line(nanmean(TBL.LocVar_post1(GIXrightIN)),'line_width',2,'style','-','color',clrs(2,:))
plot_ref_line(nanmean(TBL.LocVar_post2(GIXrightIN)),'line_width',2,'style','-','color',clrs(3,:))
legend('Baseline','L-DOPA','Ketamine');legend boxoff
pubify_figure_axis

subplot(2,1,2)
p1 = signrank(TBL.LocVar_post1(GIXleftIN) -TBL.LocVar_base(GIXleftIN));
mn1 = nanmean(TBL.LocVar_post1(GIXleftIN) -TBL.LocVar_base(GIXleftIN));
p2 = signrank(TBL.LocVar_post2(GIXleftIN) -TBL.LocVar_base(GIXleftIN));
mn2 = nanmean(TBL.LocVar_post2(GIXleftIN) -TBL.LocVar_base(GIXleftIN));
p2p1 = signrank(TBL.LocVar_post2(GIXleftIN) -TBL.LocVar_post1(GIXleftIN));
mn2m1 = nanmean(TBL.LocVar_post2(GIXleftIN) -TBL.LocVar_post1(GIXleftIN));

histogram_cowen({TBL.LocVar_base(GIXleftIN) TBL.LocVar_post1(GIXleftIN) TBL.LocVar_post2(GIXleftIN)},.1 )
title(sprintf('Left hemisphere m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f, m2m1=%1.2f, p2p1=%1.3f',mn1,p1,mn2,p2,mn2m1,p2mp1))
xlabel('Local Variance')
plot_ref_line(nanmean(TBL.LocVar_base(GIXleftIN)),'line_width',2,'style','-','color',clrs(1,:))
plot_ref_line(nanmean(TBL.LocVar_post1(GIXleftIN)),'line_width',2,'style','-','color',clrs(2,:))
plot_ref_line(nanmean(TBL.LocVar_post2(GIXleftIN)),'line_width',2,'style','-','color',clrs(3,:))
legend('Baseline','L-DOPA','Ketamine');legend boxoff
pubify_figure_axis