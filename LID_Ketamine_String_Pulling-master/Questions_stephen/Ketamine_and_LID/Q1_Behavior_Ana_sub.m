%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEHAVIOR Analysis.
% Assumes Q1_What_does_ketamine_do_firing_rate_Ana.mat is loaded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
bar(ALL.SES.Speed_Before_Early_Late)
xlabel('Recording Session')
ylabel('Movement Speed pix/sec')
pubify_figure_axis
legend(LBLS)
legend boxoff

subplot(2,1,2)
bar(ALL.SES.Jerk_Before_Early_Late)
xlabel('Recording Session')
ylabel('Jerk(|d3|)')
pubify_figure_axis

figure
Error_bars(ALL.SES.Speed_Before_Early_Late)
set(gca,'XTickLabel',LBLS)
ylabel('Pixels per second')
title('Movement before and after ketamine.')


%% ANOVA not appropriate because within subject.
figure
subplot(2,1,1)
Error_bars(ALL.SES.Speed_Before_Early_Late(:,2:3) - ALL.SES.Speed_Before_Early_Late(:,1))
set(gca,'XTickLabel',{'Post1', 'Post2'})
ylabel('Change in speed (pixels/sec)')
[~,p1] = ttest(ALL.SES.Speed_Before_Early_Late(:,2) - ALL.SES.Speed_Before_Early_Late(:,1));
[~,p2] = ttest(ALL.SES.Speed_Before_Early_Late(:,3) - ALL.SES.Speed_Before_Early_Late(:,1));
title(sprintf('Speed: p1=%0.3f,p2=%0.3f',p1,p2))

subplot(2,1,2)
Error_bars(ALL.SES.Jerk_Before_Early_Late(:,2:3) - ALL.SES.Jerk_Before_Early_Late(:,1))
set(gca,'XTickLabel',{'Post1', 'Post2'})
ylabel('Change in jerk')
[~,p11] = ttest(ALL.SES.Jerk_Before_Early_Late(:,2) - ALL.SES.Jerk_Before_Early_Late(:,1));
[~,p22] = ttest(ALL.SES.Jerk_Before_Early_Late(:,3) - ALL.SES.Jerk_Before_Early_Late(:,1));
title(sprintf('Jerk: p1=%0.3f,p2=%0.3f',p11,p22))