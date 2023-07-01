%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q2_ketamine_m1_vs_striatum_single_unit';
GP = LK_Globals;
GP.Analysis_Dir = 'C:\Temp\TempAnaResults';
PLOT_IT = true;
ses_to_ana = 'Q2_ketamine_m1_vs_striatum_single_unit';
code_dir = fileparts(which('Q2_ketamine_m1_vs_striatum_single_unit'));
adir = fullfile(GP.Analysis_Dir,ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
cm = lines(5);

%ses = {};
STR = [];
STR.group = [];
STR.day = [];
STR.animal = [];
STR.SESSIONID = [];
STR.SES.day = [];
STR.SES.animal = [];
STR.SES.group = [];

%ses = {};
M1 = [];
M1.group = [];
M1.day = [];
M1.animal = [];
M1.SESSIONID = [];
M1.SES.day = [];
M1.SES.animal = [];
M1.SES.group = [];

%%%%%%%%%%%%%%%%%%%%
STR_NRN = [];
STR_ACb = [];
STR_ACp1 = [];
STR_ACp2 = [];

% AC2b = [];
% AC2p1 = [];
% AC2p2 = [];
% 
% 
% sPSDb = [];
% sPSDp1 = [];
% sPSDp2 = [];

STRQbase = [];
STRQ = cell(3,1);


M1_NRN = [];
M1_ACb = [];
M1_ACp1 = [];
M1_ACp2 = [];

M1Qbase = [];
M1Q = cell(3,1);
ses_cnt = 1;
iEpoch =2;
nrn_cnt = 1;

%%
for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));
    [Dset.TM.BaselineUsec/60e6 diff(Dset.TM.BaselineUsec/60e6);
        Dset.TM.KetStartEndUsec/60e6 diff(Dset.TM.KetStartEndUsec/60e6);
        Dset.TM.PostInjectionUsec(1,:)/60e6 diff(Dset.TM.PostInjectionUsec(1,:)/60e6);
        Dset.TM.PostInjectionUsec(2,:)/60e6 diff(Dset.TM.PostInjectionUsec(2,:)/60e6);
        ]
    
   STR.SES.LocVar_Before_Early_Late(iF,:) = nanmean(Dset.Striatum.LocVar_Before_Early_Late);
    STRQ{1} = [STRQ{1} Dset.Striatum.Qbase];
    STRQ{2} = [STRQ{2} Dset.Striatum.Qpost1];
    STRQ{3} = [STRQ{3} Dset.Striatum.Qpost2];
    
    if PLOT_IT
        % Plot some session stuff
        figure
        bin_sec = (Dset.Striatum.Qbase_bins_uS(1,2) - Dset.Striatum.Qbase_bins_uS(1,1))/1e6;
        subplot(1,6,1:2)
        imagesc(mean(Dset.Striatum.Qbase_bins_uS,2)/60e6,[], Dset.Striatum.Qbase'/bin_sec);
        ylabel('Neuron')
        title('Baseline')
        subplot(1,6,3:4)
        imagesc(mean(Dset.Striatum.Qpost1_bins_uS,2)/60e6,[], Dset.Striatum.Qpost1'/bin_sec);
        title('Post 1')
        xlabel('min')
        subplot(1,6,5:6)
        imagesc(mean(Dset.Striatum.Qpost2_bins_uS,2)/60e6,[], Dset.Striatum.Qpost2'/bin_sec);
        title('Post 2')
        equalize_color_axes
        set(gcf,'Position',[ 41.8        503.4       1349.6        354.6])
        colorbar_label('Rate (Hz')
    end
    % Ketamine modulation
    for ii = 1:length(Dset.Striatum.FRates_base)
        STR_NRN(nrn_cnt).NeuronID = nrn_cnt;
        STR_NRN(nrn_cnt).Session = iF;
        STR_NRN(nrn_cnt).Frate_base = Dset.Striatum.FRates_base(ii);
        STR_NRN(nrn_cnt).Frate_post1 = Dset.Striatum.FRates_post1(ii);
        STR_NRN(nrn_cnt).Frate_post2 = Dset.Striatum.FRates_post2(ii);
        STR_NRN(nrn_cnt).Frate_post1mbase = Dset.Striatum.FRates_post1(ii) -  Dset.Striatum.FRates_base(ii);
        STR_NRN(nrn_cnt).Frate_post2mbase = Dset.Striatum.FRates_post2(ii) -  Dset.Striatum.FRates_base(ii);
        
        STR_NRN(nrn_cnt).LocVar_base = Dset.Striatum.LocVar_Before_Early_Late(ii,1);
        STR_NRN(nrn_cnt).LocVar_post1  = Dset.Striatum.LocVar_Before_Early_Late(ii,2);
        STR_NRN(nrn_cnt).LocVar_post2  = Dset.Striatum.LocVar_Before_Early_Late(ii,3);
        STR_NRN(nrn_cnt).LocVar_post1mbase = STR_NRN(nrn_cnt).LocVar_post1 -  STR_NRN(nrn_cnt).LocVar_base;
        STR_NRN(nrn_cnt).LocVar_post2mbase = STR_NRN(nrn_cnt).LocVar_post2 -  STR_NRN(nrn_cnt).LocVar_base;
        
        STR_ACb(nrn_cnt,:) = squeeze(Dset.Striatum.AC(ii,:,1));
        STR_ACp1(nrn_cnt,:) = squeeze(Dset.Striatum.AC(ii,:,2));
        STR_ACp2(nrn_cnt,:) = squeeze(Dset.Striatum.AC(ii,:,3));
        
        
        nrn_cnt = nrn_cnt + 1;
    end
    
end        

%%
nrn_cnt = 1;

for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));
    [Dset.TM.BaselineUsec/60e6 diff(Dset.TM.BaselineUsec/60e6);
        Dset.TM.KetStartEndUsec/60e6 diff(Dset.TM.KetStartEndUsec/60e6);
        Dset.TM.PostInjectionUsec(1,:)/60e6 diff(Dset.TM.PostInjectionUsec(1,:)/60e6);
        Dset.TM.PostInjectionUsec(2,:)/60e6 diff(Dset.TM.PostInjectionUsec(2,:)/60e6);
        ]
    
   M1.SES.LocVar_Before_Early_Late(iF,:) = nanmean(Dset.Motor.LocVar_Before_Early_Late);
    M1Q{1} = [M1Q{1} Dset.Motor.Qbase];
    M1Q{2} = [M1Q{2} Dset.Motor.Qpost1];
    M1Q{3} = [M1Q{3} Dset.Motor.Qpost2];
    
    if PLOT_IT
        % Plot some session stuff
        figure
        bin_sec = (Dset.Motor.Qbase_bins_uS(1,2) - Dset.Motor.Qbase_bins_uS(1,1))/1e6;
        subplot(1,6,1:2)
        imagesc(mean(Dset.Motor.Qbase_bins_uS,2)/60e6,[], Dset.Motor.Qbase'/bin_sec);
        ylabel('Neuron')
        title('Baseline')
        subplot(1,6,3:4)
        imagesc(mean(Dset.Motor.Qpost1_bins_uS,2)/60e6,[], Dset.Motor.Qpost1'/bin_sec);
        title('Post 1')
        xlabel('min')
        subplot(1,6,5:6)
        imagesc(mean(Dset.Motor.Qpost2_bins_uS,2)/60e6,[], Dset.Motor.Qpost2'/bin_sec);
        title('Post 2')
        equalize_color_axes
        set(gcf,'Position',[ 41.8        503.4       1349.6        354.6])
        colorbar_label('Rate (Hz')
    end
    % Ketamine modulation
    for ii = 1:length(Dset.Striatum.FRates_base)
        M1_NRN(nrn_cnt).NeuronID = nrn_cnt;
        M1_NRN(nrn_cnt).Session = iF;
        M1_NRN(nrn_cnt).Frate_base = Dset.Motor.FRates_base(ii);
        M1_NRN(nrn_cnt).Frate_post1 = Dset.Motor.FRates_post1(ii);
        M1_NRN(nrn_cnt).Frate_post2 = Dset.Motor.FRates_post2(ii);
        M1_NRN(nrn_cnt).Frate_post1mbase = Dset.Motor.FRates_post1(ii) -  Dset.Motor.FRates_base(ii);
        M1_NRN(nrn_cnt).Frate_post2mbase = Dset.Motor.FRates_post2(ii) -  Dset.Motor.FRates_base(ii);
        
        M1_NRN(nrn_cnt).LocVar_base = Dset.Motor.LocVar_Before_Early_Late(ii,1);
        M1_NRN(nrn_cnt).LocVar_post1  = Dset.Motor.LocVar_Before_Early_Late(ii,2);
        M1_NRN(nrn_cnt).LocVar_post2  = Dset.Motor.LocVar_Before_Early_Late(ii,3);
        M1_NRN(nrn_cnt).LocVar_post1mbase = M1_NRN(nrn_cnt).LocVar_post1 -  M1_NRN(nrn_cnt).LocVar_base;
        M1_NRN(nrn_cnt).LocVar_post2mbase = M1_NRN(nrn_cnt).LocVar_post2 -  M1_NRN(nrn_cnt).LocVar_base;
        
        M1_ACb(nrn_cnt,:) = squeeze(Dset.Motor.AC(ii,:,1));
        M1_ACp1(nrn_cnt,:) = squeeze(Dset.Motor.AC(ii,:,2));
        M1_ACp2(nrn_cnt,:) = squeeze(Dset.Motor.AC(ii,:,3));
        
        
        nrn_cnt = nrn_cnt + 1;
    end
    
end 

STR_TBL = struct2table(STR_NRN);
writetable(STR_TBL,'C:\Temp\Q2STR.csv')

M1_TBL = struct2table(M1_NRN);
writetable(M1_TBL,'C:\Temp\Q2M1.csv')

%%
figure
[~,six] = sort(mean(STRQ{1},1));
lbls = {'Striatum Base' 'Straitum Post 1' 'Straitum Post 2'};
for ii = 1:3
    subplot(1,3,ii)
    imagesc((1:Cols(STRQ{ii}))*bin_sec/60,[], STRQ{ii}(:,six)'/bin_sec)
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
    M = (STRQ{ii+1}-mean(STRQ{1}))'/bin_sec;
    [~,six] = sort(mean(M,2));
    subplot(1,2,ii)
    imagesc((1:Cols(STRQ{ii}))*bin_sec/60,[], M(six,:))
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
binsize_sec = diff(Dset.Striatum.Qbase_bins_uS(1,:))/1e6;
STR_mn_b = mean(STRQ{1},1)/binsize_sec;
STR_mn_p1 = mean(STRQ{2},1)/binsize_sec;
STR_mn_p2 = mean(STRQ{3},1)/binsize_sec;

figure
subplot(2,2,1)
plot(STR_mn_p1,STR_mn_b,'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(STR_mn_p2, STR_mn_b,'bo')
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
STR_p1 = signrank(STR_mn_p1-STR_mn_b);
STR_mn1 = mean(STR_mn_p1-STR_mn_b);
STR_p2 = signrank(STR_mn_p2-STR_mn_b);
STR_mn2 = mean(STR_mn_p2-STR_mn_b);

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({STR_mn_p1-STR_mn_b STR_mn_p2-STR_mn_b},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',STR_mn1,STR_p1,STR_mn2,STR_p2))

%%
%%
figure
[~,six] = sort(mean(M1Q{1},1));
lbls = {'Motor Base' 'Motor Post 1' 'Motor Post 2'};
for ii = 1:3
    subplot(1,3,ii)
    imagesc((1:Cols(M1Q{ii}))*bin_sec/60,[], M1Q{ii}(:,six)'/bin_sec)
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
    M = (M1Q{ii+1}-mean(M1Q{1}))'/bin_sec;
    [~,six] = sort(mean(M,2));
    subplot(1,2,ii)
    imagesc((1:Cols(M1Q{ii}))*bin_sec/60,[], M(six,:))
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
binsize_sec = diff(Dset.Motor.Qbase_bins_uS(1,:))/1e6;
M1_mn_b = mean(M1Q{1},1)/binsize_sec;
M1_mn_p1 = mean(M1Q{2},1)/binsize_sec;
M1_mn_p2 = mean(M1Q{3},1)/binsize_sec;

figure
subplot(2,2,1)
plot(M1_mn_p1,M1_mn_b,'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(M1_mn_p2, M1_mn_b,'bo')
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
M1_p1 = signrank(M1_mn_p1-M1_mn_b);
M1_mn1 = mean(M1_mn_p1-M1_mn_b);
M1_p2 = signrank(M1_mn_p2-M1_mn_b);
M1_mn2 = mean(M1_mn_p2-M1_mn_b);

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({M1_mn_p1-M1_mn_b M1_mn_p2-M1_mn_b},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',M1_mn1,M1_p1,M1_mn2,M1_p2))