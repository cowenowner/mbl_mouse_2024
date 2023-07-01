%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse all sessions data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q1_sal_vs_ketamine_single_unit_Ana';
GP = LK_Globals;
GP.Analysis_Dir = 'C:\Temp\TempAnaResults';
PLOT_IT = true;
ses_to_ana = 'Q1_sal_vs_ketamine_single_unit';
code_dir = fileparts(which('Q1_sal_vs_ketamine_single_unit'));
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
ACsal = [];
ACket = [];

AC2b = [];
AC2sal = [];
AC2ket = [];


sPSDb = [];
sPSDsal = [];
sPSDket = [];

AllQbase = [];
AllQ = cell(3,1);
ses_cnt = 1;
iEpoch =2;
nrn_cnt = 1;

for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));
    [Dset.TM.BaselineUsec/60e6 diff(Dset.TM.BaselineUsec/60e6);
        Dset.TM.KetStartEndUsec/60e6 diff(Dset.TM.KetStartEndUsec/60e6);
        Dset.TM.SalStartEndUsec/60e6 diff(Dset.TM.SalStartEndUsec/60e6);
        Dset.TM.SalPostInjectionUsec(1,:)/60e6 diff(Dset.TM.SalPostInjectionUsec(1,:)/60e6);
        Dset.TM.KetPostInjectionUsec(1,:)/60e6 diff(Dset.TM.KetPostInjectionUsec(1,:)/60e6);
        ]
    ALL.SES.Speed_Base_Sal_Ket(iF,:) = Dset.Speed_Base_Sal_Ket;
    ALL.SES.Jerk_Base_Sal_Ket(iF,:) = Dset.Jerk_Base_Sal_Ket;
    ALL.SES.JerkPC1_Base_Sal_Ket(iF,:) = Dset.JerkPC1_Base_Sal_Ket;
    speed_change = Dset.Speed_Base_Sal_Ket(2)-Dset.Speed_Base_Sal_Ket(1);
    jerk_change = Dset.Jerk_Base_Sal_Ket(2)-Dset.Jerk_Base_Sal_Ket(1);
    
    % the cool net level analyses.
    ALL.SES.LocVar_Base_Sal_Ket(iF,:) = nanmean(Dset.LocVar_Base_Sal_Ket);
    AllQ{1} = [AllQ{1} Dset.Qbase];
    AllQ{2} = [AllQ{2} Dset.Qsal];
    AllQ{3} = [AllQ{3} Dset.Qket];
    if PLOT_IT
        % Plot some session stuff
        figure
        bin_sec = (Dset.Qbase_bins_uS(1,2) - Dset.Qbase_bins_uS(1,1))/1e6;
        subplot(1,6,1:2)
        imagesc(mean(Dset.Qbase_bins_uS,2)/60e6,[], Dset.Qbase'/bin_sec);
        ylabel('Neuron')
        title('Baseline')
        subplot(1,6,3:4)
        imagesc(mean(Dset.Qsal_bins_uS,2)/60e6,[], Dset.Qsal'/bin_sec);
        title('Saline')
        xlabel('min')
        subplot(1,6,5:6)
        imagesc(mean(Dset.Qket_bins_uS,2)/60e6,[], Dset.Qket'/bin_sec);
        title('Ketamine')
        equalize_color_axes
        set(gcf,'Position',[ 41.8        503.4       1349.6        354.6])
        colorbar_label('Rate (Hz')
    end
    
    for ii = 1:length(Dset.FRates_base)
        NRN(nrn_cnt).NeuronID = nrn_cnt;
        NRN(nrn_cnt).Session = iF;
        NRN(nrn_cnt).Speed_change = speed_change;
        NRN(nrn_cnt).Jerk_change = jerk_change;
        NRN(nrn_cnt).Frate_base = Dset.FRates_base(ii);
        NRN(nrn_cnt).Frate_sal = Dset.FRates_sal(ii);
        NRN(nrn_cnt).Frate_ket = Dset.FRates_ket(ii);
        NRN(nrn_cnt).Frate_salmbase = Dset.FRates_sal(ii) -  Dset.FRates_base(ii);
        NRN(nrn_cnt).Frate_ketmbase = Dset.FRates_ket(ii) -  Dset.FRates_base(ii);
        
        NRN(nrn_cnt).LocVar_base = Dset.LocVar_Base_Sal_Ket(ii,1);
        NRN(nrn_cnt).LocVar_sal  = Dset.LocVar_Base_Sal_Ket(ii,2);
        NRN(nrn_cnt).LocVar_ket  = Dset.LocVar_Base_Sal_Ket(ii,3);
        NRN(nrn_cnt).LocVar_salmbase = NRN(nrn_cnt).LocVar_sal -  NRN(nrn_cnt).LocVar_base;
        NRN(nrn_cnt).LocVar_ketmbase = NRN(nrn_cnt).LocVar_ket -  NRN(nrn_cnt).LocVar_base;
        
        ACb(nrn_cnt,:) = squeeze(Dset.AC(ii,:,1));
        ACsal(nrn_cnt,:) = squeeze(Dset.AC(ii,:,2));
        ACket(nrn_cnt,:) = squeeze(Dset.AC(ii,:,3));
                
        nrn_cnt = nrn_cnt + 1;
    end
    
end

% Save it to a nice table that is easy to read in R.
TBL = struct2table(NRN);
writetable(TBL,'C:\Temp\Q1SalKet.csv')

LBLS = {'base','post saline', 'post ketamine'};

figure
subplot(2,1,1)
bar(ALL.SES.Speed_Base_Sal_Ket)
xlabel('Recording Session')
ylabel('Movement Speed pix/sec')
pubify_figure_axis
legend(LBLS)
legend boxoff

subplot(2,1,2)
bar(ALL.SES.Jerk_Base_Sal_Ket)
xlabel('Recording Session')
ylabel('Jerk(|d3|)')
pubify_figure_axis

%%
figure
[~,six] = sort(mean(AllQ{1},1));
lbls = {'Base' 'Saline' 'Ketamine'};
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
mn_sal = mean(AllQ{2},1)/binsize_sec;
mn_ket = mean(AllQ{3},1)/binsize_sec;

figure
subplot(2,1,1)
plot(mn_sal,mn_b,'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Saline (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,1,2)
plot(mn_ket, mn_b,'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Ketamine (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch')

%% Autocorr analysis...
BASEIX = Dset.AC_x_ms > 120;
XIX = Dset.AC_x_ms < 60;
if PLOT_IT
    
    % Plot one at a time and PAUSE - press space to continue 
    for ii = 1:Rows(ACb)
        figure
        plot(Dset.AC_x_ms(XIX),ACb(ii,XIX),'Color',cm(1,:),'LineWidth',4)
        hold on
        plot(Dset.AC_x_ms(XIX),ACsal(ii,XIX),'Color',cm(2,:),'LineWidth',4)
        plot(Dset.AC_x_ms(XIX),ACket(ii,XIX),'Color',cm(3,:),'LineWidth',4)
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
imagesc(Dset.AC_x_ms(XIX),[],ACsal(six,XIX) - mean(ACb(six,BASEIX),2))
pubify_figure_axis
title('Saline')
subplot(1,3,3)
imagesc(Dset.AC_x_ms(XIX),[],ACket(six,XIX) - mean(ACb(six,BASEIX),2))
pubify_figure_axis
title('Ketamine')

equalize_color_axes
colorbar_label

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot_confidence_intervals(Dset.AC_x_ms(XIX),ACb(:,XIX) - mean(ACb(:,BASEIX),2),[],cm(1,:))
hold on
plot_confidence_intervals(Dset.AC_x_ms(XIX),ACsal(:,XIX) - mean(ACb(:,BASEIX),2),[],cm(2,:))
plot_confidence_intervals(Dset.AC_x_ms(XIX),ACket(:,XIX) - mean(ACb(:,BASEIX),2),[],cm(3,:))
legend_color_text(LBLS,cm)
pubify_figure_axis
xlabel('ms')
ylabel('Normalized Rate')
title('Mean Auto Corrs')

