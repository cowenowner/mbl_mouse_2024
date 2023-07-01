%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q1_Do_neurons_correlate_with_string_pulling_Ana';
GP = LK_Globals;
GP.Analysis_Dir = 'E:\Temp\TempAnaResults';
PLOT_IT = true;
ses_to_ana = 'Q1_Do_neurons_correlate_with_string_pulling';
code_dir = fileparts(which('Q1_Do_neurons_correlate_with_string_pulling'));
adir = fullfile(GP.Analysis_Dir,ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
cm = lines(5);

ses = {};
ALL = [];
ALL.group = [];
ALL.day = [];
ALL.animal = [];
ALL.SESSIONID = [];
ALL.SES.day = [];
ALL.SES.animal = [];
ALL.SES.group = [];

%%%%%%%%%%%%%%%%%%%%
NRN = [];
ACb = [];
ACp1 = [];
ACp2 = [];

AC2b = [];
AC2p1 = [];
AC2p2 = [];


sPSDb = [];
sPSDp1 = [];
sPSDp2 = [];

AllQbase = [];
AllQ = cell(3,1);
ses_cnt = 1;
iEpoch =2;
nrn_cnt = 1;
for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));
    [Dset.TM.BaselineUsec/60e6 diff(Dset.TM.BaselineUsec/60e6);
        Dset.TM.KetStartEndUsec/60e6 diff(Dset.TM.KetStartEndUsec/60e6);
        Dset.TM.PostInjectionUsec(1,:)/60e6 diff(Dset.TM.PostInjectionUsec(1,:)/60e6);
        Dset.TM.PostInjectionUsec(2,:)/60e6 diff(Dset.TM.PostInjectionUsec(2,:)/60e6);
        ]
    ALL.SES.Speed_Before_Early_Late(iF,:) = Dset.Speed_Before_Early_Late;
    ALL.SES.Jerk_Before_Early_Late(iF,:) = Dset.Jerk_Before_Early_Late;
    ALL.SES.JerkPC1_Before_Early_Late(iF,:) = Dset.JerkPC1_Before_Early_Late;
    speed_change = Dset.Speed_Before_Early_Late(2)-Dset.Speed_Before_Early_Late(1);
    jerk_change = Dset.Jerk_Before_Early_Late(2)-Dset.Jerk_Before_Early_Late(1);
    
    % the cool net level analyses.
    ALL.SES.LocVar_Before_Early_Late(iF,:) = nanmean(Dset.LocVar_Before_Early_Late);
    AllQ{1} = [AllQ{1} Dset.Qbase];
    AllQ{2} = [AllQ{2} Dset.Qpost1];
    AllQ{3} = [AllQ{3} Dset.Qpost2];
    if PLOT_IT
        % Plot some session stuff
        figure
        bin_sec = (Dset.Qbase_bins_uS(1,2) - Dset.Qbase_bins_uS(1,1))/1e6;
        subplot(1,6,1:2)
        imagesc(mean(Dset.Qbase_bins_uS,2)/60e6,[], Dset.Qbase'/bin_sec);
        ylabel('Neuron')
        title('Baseline')
        subplot(1,6,3:4)
        imagesc(mean(Dset.Qpost1_bins_uS,2)/60e6,[], Dset.Qpost1'/bin_sec);
        title('Post 1')
        xlabel('min')
        subplot(1,6,5:6)
        imagesc(mean(Dset.Qpost2_bins_uS,2)/60e6,[], Dset.Qpost2'/bin_sec);
        title('Post 2')
        equalize_color_axes
        set(gcf,'Position',[ 41.8        503.4       1349.6        354.6])
        colorbar_label('Rate (Hz')
    end
    % Ketamine modulation
    for ii = 1:length(Dset.FRates_base)
        NRN(nrn_cnt).NeuronID = nrn_cnt;
        NRN(nrn_cnt).Session = iF;
        NRN(nrn_cnt).Speed_change = speed_change;
        NRN(nrn_cnt).Jerk_change = jerk_change;
        NRN(nrn_cnt).Frate_base = Dset.FRates_base(ii);
        NRN(nrn_cnt).Frate_post1 = Dset.FRates_post1(ii);
        NRN(nrn_cnt).Frate_post2 = Dset.FRates_post2(ii);
        NRN(nrn_cnt).Frate_post1mbase = Dset.FRates_post1(ii) -  Dset.FRates_base(ii);
        NRN(nrn_cnt).Frate_post2mbase = Dset.FRates_post2(ii) -  Dset.FRates_base(ii);
        
        NRN(nrn_cnt).LocVar_base = Dset.LocVar_Before_Early_Late(ii,1);
        NRN(nrn_cnt).LocVar_post1  = Dset.LocVar_Before_Early_Late(ii,2);
        NRN(nrn_cnt).LocVar_post2  = Dset.LocVar_Before_Early_Late(ii,3);
        NRN(nrn_cnt).LocVar_post1mbase = NRN(nrn_cnt).LocVar_post1 -  NRN(nrn_cnt).LocVar_base;
        NRN(nrn_cnt).LocVar_post2mbase = NRN(nrn_cnt).LocVar_post2 -  NRN(nrn_cnt).LocVar_base;
        
        ACb(nrn_cnt,:) = squeeze(Dset.AC(ii,:,1));
        ACp1(nrn_cnt,:) = squeeze(Dset.AC(ii,:,2));
        ACp2(nrn_cnt,:) = squeeze(Dset.AC(ii,:,3));
        
        %         AC2b(nrn_cnt,:) = squeeze(Dset.AcorrSmth_Before_Early_Late(ii,:,1));
        %         AC2p1(nrn_cnt,:) = squeeze(Dset.AcorrSmth_Before_Early_Late(ii,:,2));
        %         AC2p2(nrn_cnt,:) = squeeze(Dset.AcorrSmth_Before_Early_Late(ii,:,3));
        %
        %         sPSDb(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_Late(ii,:,1));
        %         sPSDp1(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_Late(ii,:,2));
        %         sPSDp2(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_Late(ii,:,3));
        %
        %         sPSDbZ(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_LateZ(ii,:,1));
        %         sPSDp1Z(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_LateZ(ii,:,2));
        %         sPSDp2Z(nrn_cnt,:) = squeeze(Dset.SpikePSD_Before_Early_LateZ(ii,:,3));
        
        nrn_cnt = nrn_cnt + 1;
    end
    
end

% Save it to a nice table that is easy to read in R.
TBL = struct2table(NRN);
writetable(TBL,'C:\Temp\Q1Ket.csv')

LBLS = {'base','post 1', 'post 2'};

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
%%
figure
[~,six] = sort(mean(AllQ{1},1));
lbls = {'Base' 'Post 1' 'Post 2'};
for ii = 1:3
    subplot(1,3,ii)
    imagesc((1:Cols(AllQ{ii}))*bin_sec/60,[], AllQ{ii}(:,six)'/bin_sec)
    xlabel('Min')
    if ii == 1
        ylabel('Neuron')
    end
    title(lbls{ii})

    pubify_figure_axis
end
equalize_color_axes
colorbar_label

figure
for ii = 1:2
    M = (AllQ{ii+1}-mean(AllQ{1}))'/bin_sec;
    [~,six] = sort(mean(M,2));
    subplot(1,2,ii)
    imagesc((1:Cols(AllQ{ii}))*bin_sec/60,[], M(six,:))
    xlabel('Min')
    ylabel('Neuron')
    title(sprintf('Diff from base Post %d',ii))
    pubify_figure_axis
    caxis([-3 3])
    colormap(hotcold_colormap)
end
equalize_color_axes
colorbar_label
%% Mean firing rate analyses...
binsize_sec = diff(Dset.Qbase_bins_uS(1,:))/1e6;
mn_b = mean(AllQ{1},1)/binsize_sec;
mn_p1 = mean(AllQ{2},1)/binsize_sec;
mn_p2 = mean(AllQ{3},1)/binsize_sec;

figure
subplot(2,2,1)
plot(mn_p1,mn_b,'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(mn_p2, mn_b,'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch')

subplot(2,2,3:4)
p1 = signrank(mn_p1-mn_b);
mn1 = mean(mn_p1-mn_b);
p2 = signrank(mn_p2-mn_b);
mn2 = mean(mn_p2-mn_b);

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1-mn_b mn_p2-mn_b},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))
%% Popluation distance comparison (as in that Mohgaddahm paper we went over...

figure
subplot(1,2,1)
p1 = pdist([mn_b;mn_p1],'euclidean');
p2 = pdist([mn_b;mn_p2],'euclidean');
bar([p1 p2])
set(gca,'YLim',[12 19])
ylabel('Euclidean Distance from Baseline')
set(gca,'XTickLabel',{'P1-B','P2-B'})
pubify_figure_axis

subplot(1,2,2)
p1 = pdist([mn_b;mn_p1],'spearman');
p2 = pdist([mn_b;mn_p2],'spearman');
bar([p1 p2])
ylabel('Spearman r Distance from Baseline')
set(gca,'XTickLabel',{'P1-B','P2-B'})
pubify_figure_axis

%% Autocorr analysis...
BASEIX = Dset.AC_x_ms > 120;
XIX = Dset.AC_x_ms < 60;
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
[~,six] = sort(nanmax(ACb,[],2));

figure
subplot(1,3,1)
imagesc(Dset.AC_x_ms(XIX),[],ACb(six,XIX) - mean(ACb(six,BASEIX),2))
xlabel('ms')
ylabel('Neuron (sorted by peak)')
title('Baseline Norm Autocorr')
pubify_figure_axis
subplot(1,3,2)
imagesc(Dset.AC_x_ms(XIX),[],ACp1(six,XIX) - mean(ACb(six,BASEIX),2))
pubify_figure_axis
title('Post 1')
subplot(1,3,3)
imagesc(Dset.AC_x_ms(XIX),[],ACp2(six,XIX) - mean(ACb(six,BASEIX),2))
pubify_figure_axis
title('Post 2')

equalize_color_axes
colorbar_label

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot_confidence_intervals(Dset.AC_x_ms(XIX),ACb(:,XIX) - mean(ACb(:,BASEIX),2),[],cm(1,:))
hold on
plot_confidence_intervals(Dset.AC_x_ms(XIX),ACp1(:,XIX) - mean(ACb(:,BASEIX),2),[],cm(2,:))
plot_confidence_intervals(Dset.AC_x_ms(XIX),ACp2(:,XIX) - mean(ACb(:,BASEIX),2),[],cm(3,:))
legend_color_text(LBLS,cm)
pubify_figure_axis
xlabel('ms')
ylabel('Normalized Rate')
title('Mean Auto Corrs')

% XIX = (Dset.AcorrSmth_x_ms) > 0 & (Dset.AcorrSmth_x_ms) < 80 ;
% 
% figure
% [~,six] = sort(nanmax(AC2b,[],2));
% subplot(1,3,1)
% imagesc(Dset.AcorrSmth_x_ms(XIX),[],AC2b(six,XIX))
% xlabel('ms')
% ylabel('Neuron')
% pubify_figure_axis
% subplot(1,3,2)
% imagesc(Dset.AcorrSmth_x_ms(XIX),[],AC2p1(six,XIX))
% pubify_figure_axis
% subplot(1,3,3)
% imagesc(Dset.AcorrSmth_x_ms(XIX),[],AC2p2(six,XIX))
% pubify_figure_axis
% 
% equalize_color_axes
% colorbar_label
% 
% 
% figure
% plot_confidence_intervals(Dset.AcorrSmth_x_ms(XIX),AC2b(:,XIX),[],cm(1,:))
% hold on
% plot_confidence_intervals(Dset.AcorrSmth_x_ms(XIX),AC2p1(:,XIX),[],cm(2,:))
% plot_confidence_intervals(Dset.AcorrSmth_x_ms(XIX),AC2p2(:,XIX),[],cm(3,:))
% legend_color_text(LBLS,cm)
% pubify_figure_axis
% xlabel('ms')
% ylabel('r')


% 
% %% PSDs
% figure
% subplot(1,3,1)
% imagesc(Dset.SpikePSD_fqs,[],sPSDb)
% xlabel('Hz')
% ylabel('Neuron')
% pubify_figure_axis
% 
% subplot(1,3,2)
% imagesc(Dset.SpikePSD_fqs,[],sPSDp1)
% xlabel('Hz')
% ylabel('Neuron')
% pubify_figure_axis
% 
% subplot(1,3,3)
% imagesc(Dset.SpikePSD_fqs,[],sPSDp2)
% xlabel('Hz')
% ylabel('Neuron')
% pubify_figure_axis
% 
% 
% %% PSDsZ
% figure
% subplot(2,3,1)
% imagesc(Dset.SpikePSD_fqs,[],sPSDbZ)
% caxis([-2 5])
% xlabel('Hz')
% ylabel('Neuron')
% pubify_figure_axis
% 
% subplot(2,3,2)
% imagesc(Dset.SpikePSD_fqs,[],sPSDp1Z)
% xlabel('Hz')
% ylabel('Neuron')
% pubify_figure_axis
% caxis([-2 5])
% 
% subplot(2,3,3)
% imagesc(Dset.SpikePSD_fqs,[],sPSDp2Z)
% xlabel('Hz')
% ylabel('Neuron')
% pubify_figure_axis
% caxis([-2 5])
% colorbar_label('Z')
% 
% subplot(2,3,4)
% plot_confidence_intervals(Dset.SpikePSD_fqs,sPSDbZ);
% set(gca,'YLim',[-1,3])
% plot_horiz_line_at_zero
% ylabel('std')
% subplot(2,3,5)
% plot_confidence_intervals(Dset.SpikePSD_fqs,sPSDp1Z);
% set(gca,'YLim',[-1,3])
% plot_horiz_line_at_zero
% xlabel('Hz')
% subplot(2,3,6)
% plot_confidence_intervals(Dset.SpikePSD_fqs,sPSDp2Z);
% set(gca,'YLim',[-1,3])
% plot_horiz_line_at_zero
