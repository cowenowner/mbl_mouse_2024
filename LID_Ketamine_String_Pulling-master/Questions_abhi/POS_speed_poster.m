%% loop for all sessions
data_dir = 'E:\SfN_LFP';
sFreq = 50;
baseline_min = [-60 -40];
ana_range_min = [-60 60];
ana_t_pos_sec = ana_range_min(1)*60:1/sFreq:ana_range_min(2)*60;
ALL.Speed = [];
conditions_to_compare = {'POS_speed_Sal_Ket_LID' 'POS_speed_Sal_Ket_control' };

for iD = 1:length(conditions_to_compare)
    files = dir(fullfile(data_dir,conditions_to_compare{iD},'*.mat'));
    for iF = 1:length(files)
        D = load(fullfile(data_dir,conditions_to_compare{iD},files(iF).name));
        D.POS.speed_pix_sec = Speed_from_xy([D.POS.Time_uS/1e6 double(D.POS.Red_xy)]);
        D.t_sec = D.POS.Time_uS/1e6;
        if any(contains(fieldnames(D.POS),'Sal_Start_min'))
            D.t_sec = D.t_sec - D.POS.Sal_Start_min*60;
        else
            D.t_sec = D.t_sec - D.POS.Ket_Start_min*60;
        end
        POS = single(interp1(D.t_sec,D.POS.speed_pix_sec,ana_t_pos_sec));
        POS = [ana_t_pos_sec(:) POS(:)];
        IX = POS(:,1) > baseline_min(1)*60 & POS(:,1) < baseline_min(2)*60;
        mn = nanmean(POS(IX,2)');
        sd = nanstd(POS(IX,2)');
        POS(:,3) = (POS(:,2)-mn)./sd;
        ALL.Speed{iD}(iF,:) = POS(:,3)';
    end
end

%% doing stats
ket_eff_win = [2 30];
IX_ket = POS(:,1) > ket_eff_win(1)*60 & POS(:,1) < ket_eff_win(2)*60;

for ii = 1:length(conditions_to_compare)
    PS_allses = squeeze(ALL.Speed{ii}(:,IX_ket));
    PS_mn_allses = nanmean(PS_allses');
    PS_allcondition{ii} = PS_mn_allses';
end

[h,pval] = ranksum(PS_allcondition{1},PS_allcondition{2});

% for plotting
for ii = 1:length(conditions_to_compare)
    PS_allses = squeeze(ALL.Speed{ii}(:,IX_ket));
    PS_mn_allses = nanmean(PS_allses');
    PS_allcondition{ii} = PS_mn_allses;
end

PS_concat = cell2mat(PS_allcondition);
g1 = 1*ones(1,8);
g2 = 2*ones(1,10);
group = [g1,g2];

figure; boxplot(PS_concat',group')
% plotting

X1=1;
X2=2;
figure;

boxplot(PS_concat',group','colors','r','widths',0.25)

axis tight

figure
violin([PS_concat(1:8) PS_concat(9:end)])

%% movement 
ana_range_min = [100 310];

IMU = LK_Load_and_Process_IMU('Inertial_data.mat');
POS = LK_Load_and_Clean_POS('POS.mat');
ana_t_pos_sec = ana_range_min(1)*60:1/IMU.sFreq:ana_range_min(2)*60;
IX = POS.Time_uS > ana_range_min(1)*60e6 & POS.Time_uS < ana_range_min(2)*60e6;

s = movmedian(POS.speed,3000,'omitnan');
s = conv_filter(POS.speed,hanning(1000));
figure
plot(POS.Time_uS/60e6,s)
axis tight
figure;plot(IX)

% restrict speed to ana range
t = POS.Time_uS(IX,:);
s_range = s(IX,:);
figure
plot(t/60e6,s_range)