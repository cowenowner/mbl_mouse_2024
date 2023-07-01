function [OUT] = LK_Find_best_LFP_for_a_band_Abhi()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine which channel or channel combo has the strongest
% oscillation
%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global DIRS
% if nargin < 1
%     sig_name = 'gamma_80';
% end
OUT = [];
% [GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();
SES = LK_Session_Info();
OUT.SES = SES;
OUT.aborted = false;
fqs = 1:.5:170;
LFP_analysis_dir = fullfile('H:\LFP_analysis',SES.rat_str,SES.session_str);
if ~exist(LFP_analysis_dir,'dir')
    mkdir(LFP_analysis_dir)
end
% get the LFP directory
DATA_DIR = pwd;
tmp = dir(fullfile(DATA_DIR,'Recording*'));
lfp_root_dir = fullfile(DATA_DIR,tmp(1).name);
lfp_dir = fullfile(DATA_DIR,tmp(1).name,'LFP');

% lfp_dir = fullfile{LFP_Dir,SES.rat_str,SES.session_str,'LFP');
%lfp_dir = fullfile(DIRS.LFP_Dir,SES.rat_str,SES.session_str,'LFP');
[~,LFP_files] = find_files(fullfile(lfp_dir,'amp*.mat'));

% Seperate the right and left hemispheres
% Get the TT and channel number for LFP files
ch_num = [];
for ii = 1:length(LFP_files)
    Part1 = split(LFP_files{ii},'_');
    Part2 = split(Part1{1},'-');
    ch_num(ii) = str2double(Part2{3});
end
% Load LFP to tetrode conversion excel
LtT = readtable('../LFP_to_Tetrode.xlsx');
TT_num = [];
for ii = 1:length(ch_num)
    IX = LtT.LFPChannel == ch_num(ii);
    if any(IX)
        TT_num(ii) = LtT.Tetrode(IX);
    end
end

% Get the depths of the TT
load(fullfile(lfp_root_dir,'Meta_data.mat'))
[META.recdatestr] = LK_Determine_Recording_Time_From_Datadir;

DEPTHS = LK_Load_Depths('..',META.recdatestr);
Depth_lf = [];
for ii = 1:length(TT_num)
    IX = DEPTHS(:,1) == TT_num(ii);
    if any(IX)
        Depth_lf(ii) = DEPTHS(IX,2);
    end
end

T2C = readtable('../Tetrode_2d_Cordinates.xlsx');
Location_lf = [];
for ii = 1:length(TT_num)
    IX = T2C.Tetrode == TT_num(ii);
    if any(IX)
        Location_lf(ii) = T2C.MLmm(IX);
    end
end

LFP_files_info = [ch_num' TT_num' Depth_lf' Location_lf'];

OUT.LFP_ch_TT_depth_info = LFP_files_info;


E = LK_Load_Events_From_Excel('Event_times.xlsx');
% switch sig_name
% %     case {'gamma_80' 'gamma_50' 'beta'}
%     case 'gamma_80'
        if any(E.EventID=='LDOPAInjectionStart')
            t_start_end_min(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')+38;
            t_start_end_min(2) = t_start_end_min(1) + 20;
            sig_name = 'gamma_80';
        elseif  any(E.EventID=='KetInjectionStart')
            t_start_end_min(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')+2;
            t_start_end_min(2) = t_start_end_min(1) + 20;
            sig_name = 'gamma_50';
        end
        
%     case 'beta'
%         if any(E.EventID=='TreadmillStart')
%             disp('Found Treadmill')
%             tmp1 = E.MinFromStart(E.EventID == 'TreadmillStart');
%             tmp2 = E.MinFromStart(E.EventID == 'TreadmillStop');
%             
%             t_start_end_min(1) = tmp1(1);
%             t_start_end_min(2) = tmp2(1);
%             
%                         t_start_end_min(1) = tmp1(1)-20;
%                         t_start_end_min(2) = tmp1(1)-2;
%             
%         else
%         end
       
%     otherwise
%         if any(E.EventID=='KetInjectionStart')
%             disp('Found Ketamine')
%             
%             t_start_end_min(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')+10;
%             t_start_end_min(2) = t_start_end_min(1) + 20;
%         end
        
% end
OUT.sig_name = sig_name;
t_start_end_uS = t_start_end_min*60e6;


for iF = 1:length(LFP_files)
    LFP = LK_Load_and_Clean_LFP(lfp_dir, LFP_files{iF});
    if iF == 1
        IX = LFP.t_uS > t_start_end_uS(1) & LFP.t_uS < t_start_end_uS(2);
        L = zeros(sum(IX),length(LFP_files));
        L(:,1) = LFP.LFP(IX);
    else
        L(:,iF) = LFP.LFP(IX);
    end
end
switch sig_name
    case 'gamma_50'
        filts = SPEC_create_filters({'gamma_50' 'beta2'},LFP.sFreq);
    case 'beta'
        filts = SPEC_create_filters({[25 37] [15 21]},LFP.sFreq);
    case 'gamma_80'
        filts = SPEC_create_filters({'gamma_80' [65 75]},LFP.sFreq);
    case 'HFO'
        filts = SPEC_create_filters({'HFO' [110 125]},LFP.sFreq);
end

OUT.filts = filts;

% Get Right hem and Left hem M1 channels
IX_R = ch_num > 31 & Depth_lf < 2000;

IX_L = ch_num < 32 & Depth_lf < 2000;

IX_RS = ch_num > 31 & Depth_lf > 2800;
IX_LS = ch_num < 32 & Depth_lf > 2800; 

[OUT.INFO_R_M1] = SPEC_find_best_channel_combo(L(:,IX_R),filts{1},filts{2});
[OUT.INFO_L_M1] = SPEC_find_best_channel_combo(L(:,IX_L),filts{1},filts{2});

% best combo for right hemisphere M1
LFP_files_R = LFP_files(IX_R);
OUT.best_non_reref_Right_M1 = LFP_files_R{OUT.INFO_R_M1.best_non_reref};
OUT.best_reref_combo_Right_M1{1} = LFP_files_R{OUT.INFO_R_M1.best_reref_combo(1)};
OUT.best_reref_combo_Right_M1{2} = LFP_files_R{OUT.INFO_R_M1.best_reref_combo(2)};
OUT.best_non_reref_ratio_Right_M1 = LFP_files_R{OUT.INFO_R_M1.best_non_reref_ratio};
OUT.best_non_reref_prop_Right_M1 = LFP_files_R{OUT.INFO_R_M1.best_non_reref_prop};


LFP_info_R = LFP_files_info(IX_R,:);  
dist_from_best_R = LFP_info_R(OUT.INFO_R_M1.best_non_reref_ratio,4) - LFP_info_R(:,4);
IX_dist_700_R = abs(dist_from_best_R) > .7;
if sum(IX_dist_700_R) > 1
    [max_dist, max_id] = max(abs(dist_from_best_R(IX_dist_700_R)));
    OUT.best_reref_700uM_apart_Right_M1 = LFP_files_R{max_id};
    OUT.best_reref_dist_apart_Right_M1 = max_dist;
elseif sum(IX_dist_700_R) == 1
    OUT.best_reref_700uM_apart_Right_M1 = LFP_files_R{IX_dist_700_R};
    OUT.best_reref_dist_apart_Right_M1 = dist_from_best_R(IX_dist_700_R);
end
    
% best combo for left hemisphere M1
LFP_files_L = LFP_files(IX_L);
OUT.best_non_reref_Left_M1 = LFP_files_L{OUT.INFO_L_M1.best_non_reref};
OUT.best_reref_combo_Left_M1{1} = LFP_files_L{OUT.INFO_L_M1.best_reref_combo(1)};
OUT.best_reref_combo_Left_M1{2} = LFP_files_L{OUT.INFO_L_M1.best_reref_combo(2)};
OUT.best_non_reref_ratio_Left_M1 = LFP_files_L{OUT.INFO_L_M1.best_non_reref_ratio};
OUT.best_non_reref_prop_Left_M1 = LFP_files_L{OUT.INFO_L_M1.best_non_reref_prop};

LFP_info_L = LFP_files_info(IX_L,:);  
dist_from_best_L = LFP_info_L(OUT.INFO_L_M1.best_non_reref_ratio,4) - LFP_info_L(:,4);
IX_dist_700_L = abs(dist_from_best_L) > .7;
if sum(IX_dist_700_L) > 1
    [max_dist, max_id] = max(abs(dist_from_best_L(IX_dist_700_L)));
    OUT.best_reref_700uM_apart_Left_M1 = LFP_files_L{max_id};
    OUT.best_reref_dist_apart_Left_M1 = max_dist;
elseif sum(IX_dist_700_L) == 1
    OUT.best_reref_700uM_apart_Left_M1 = LFP_files_L{IX_dist_700_L};
    OUT.best_reref_dist_apart_Left_M1 = dist_from_best_L(IX_dist_700_L);
end

% Get the Striatum stuff
if any(IX_RS)
    [OUT.INFO_R_STR] = SPEC_find_best_channel_combo(L(:,IX_RS),filts{1},filts{2});    
    LFP_files_RS = LFP_files(IX_RS);
    OUT.best_non_reref_Right_STR = LFP_files_RS{OUT.INFO_R_STR.best_non_reref};
    OUT.best_non_reref_ratio_Right_STR = LFP_files_RS{OUT.INFO_R_STR.best_non_reref_ratio};
    OUT.best_non_reref_prop_Right_STR = LFP_files_RS{OUT.INFO_R_STR.best_non_reref_prop};

    copyfile(fullfile(lfp_dir,OUT.best_non_reref_ratio_Right_STR),fullfile(LFP_analysis_dir,['Right_STR_best_nonreref_ratio' OUT.best_non_reref_ratio_Right_STR '.mat']))
else
end

if any(IX_LS)
    [OUT.INFO_L_STR] = SPEC_find_best_channel_combo(L(:,IX_LS),filts{1},filts{2});
    LFP_files_LS = LFP_files(IX_LS);
    OUT.best_non_reref_Left_STR = LFP_files_LS{OUT.INFO_L_STR.best_non_reref};
    OUT.best_non_reref_ratio_Left_STR = LFP_files_LS{OUT.INFO_L_STR.best_non_reref_ratio};
    OUT.best_non_reref_prop_Left_STR = LFP_files_LS{OUT.INFO_L_STR.best_non_reref_prop};
    copyfile(fullfile(lfp_dir,OUT.best_non_reref_ratio_Left_STR),fullfile(LFP_analysis_dir,['Left_STR_best_nonreref_ratio' OUT.best_non_reref_ratio_Left_STR '.mat']))
else
end

% make the last col the re-reffed data.
% RR = L(:,OUT.INFO.best_reref_combo(1)) - L(:,OUT.INFO.best_reref_combo(2)) ;

% OUT
save(['best_channels_' sig_name '.mat'],'OUT')


copyfile(fullfile(lfp_dir,OUT.best_non_reref_ratio_Right_M1),fullfile(LFP_analysis_dir,['Right_M1_best_nonreref_ratio' OUT.best_non_reref_ratio_Right_M1 '.mat']))
copyfile(fullfile(lfp_dir,OUT.best_non_reref_ratio_Left_M1),fullfile(LFP_analysis_dir,['Left_M1_best_nonreref_ratio' OUT.best_non_reref_ratio_Left_M1 '.mat']))
if isfield(OUT,'best_reref_700uM_apart_Right_M1')
    copyfile(fullfile(lfp_dir,OUT.best_reref_700uM_apart_Right_M1),fullfile(LFP_analysis_dir,['Right_M1_700uM_apart_' OUT.best_reref_700uM_apart_Right_M1 '.mat']))
else
end
if isfield(OUT,'best_reref_700uM_apart_Left_M1')
    copyfile(fullfile(lfp_dir,OUT.best_reref_700uM_apart_Left_M1),fullfile(LFP_analysis_dir,['Left_M1_700uM_apart_' OUT.best_reref_700uM_apart_Left_M1 '.mat']))
else
end

copyfile(fullfile(DATA_DIR,'*.mat'),LFP_analysis_dir);
copyfile(fullfile(DATA_DIR,'*.xlsx'),LFP_analysis_dir);
copyfile(fullfile(DATA_DIR,'*.docx'),LFP_analysis_dir);
copyfile(fullfile(lfp_root_dir,'*.mat'),LFP_analysis_dir);
copyfile(fullfile(lfp_root_dir,'*.rhd'),LFP_analysis_dir);
copyfile(fullfile(lfp_root_dir,'*.pos'),LFP_analysis_dir);



%%
% figure
% subplot(2,1,1)
% pwelch(L(:,OUT.INFO.best_reref_combo(1)),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
% grid off
% pubify_figure_axis
% title('')
% ylabel('')
% subplot(2,1,2)
% pwelch(L(:,OUT.INFO.best_reref_combo(2)),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
% grid off
% title('')
% ylabel('')
% pubify_figure_axis
% equalize_y_axes
% 
% figure
% pwelch(RR,LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
% title('re reffed signal')
% 
% plot_me(OUT,LFP_files,sig_name)
% plot_all_LFP(L,LFP,LFP_files,sig_name)
% % plot_all_LFP_spec(LFP,LFP_files,sig_name) - needs cleaning
% % Now that this was selected, show the spectrogram.
% files = {OUT.best_non_reref OUT.best_reref_combo{1} OUT.best_reref_combo{2} 'reref'};
% chnls = [OUT.INFO.best_non_reref OUT.INFO.best_reref_combo Cols(L)];
% %%
% figure
% for iF = 1:length(files)
%     if strcmp(files{iF},'reref')
%         a = LK_Load_and_Clean_LFP(lfp_dir, OUT.best_reref_combo{1});
%         LFP = LK_Load_and_Clean_LFP(lfp_dir, OUT.best_reref_combo{2});
%         LFP.LFP = a.LFP-LFP.LFP;
%     else
%         LFP = LK_Load_and_Clean_LFP(lfp_dir, files{iF});
%         
%     end
%     % whiten: Does not work so great.
%     %     LFP.LFP = SPEC_whiten_data_AR(LFP.LFP,3);
%     subplot(length(files),1,iF)
%     %      spectrogram(LFP.LFP,LFP.sFreq*5,LFP.sFreq*3,fqs, LFP.sFreq,'yaxis');
%     spec_win_sec = 5;
%     [S,w,t] = spectrogram(LFP.LFP,LFP.sFreq*spec_win_sec,LFP.sFreq*(spec_win_sec/2),fqs, LFP.sFreq);
%     sFreq_of_spectrogram = 1/median(diff(t));
%     t_min = t/60;
%     S = 10*log10(abs(S))';
%     %     p = prctile(S,[.5 99.5]);
%     %     BIX = S<p(1,:);
%     %     BIX = BIX | S>p(2,:);
%     %     S(BIX) = nan;
% %     Sn = movmedian(S,30,'omitnan');
%     han_size_sec = 40;
%     han_size_points = han_size_sec*sFreq_of_spectrogram;
%     Sn = conv_filter(S,hanning(han_size_points)/sum(hanning(han_size_points)));
%     
%     BASEIX = t_min < 10;
%     Sb = Sn-nanmean(Sn(BASEIX,:),1);
%     Sb = Sb./nanstd(Sn(BASEIX,:),[],1);
% 
%     imagesc(t_min,fqs,Sn');
%     
%     axis xy
%     mkrs = cat2cell(E.EventID);
%     plot_markers(E.MinFromStart,1)
%     title(files{iF})
%     ylabel('Hz')
%     colorbar
%     caxis(prctile(Sn(:),[1 98.6]))
%     
%     if iF == length(files)
%         h = plot_markers(E.MinFromStart,1,mkrs);
%         
%         set(h,'Rotation',45)
%         set(h,'FontSize',5)
%         xlabel('min')
%     end
% end
% % POS = LK_Load_and_Clean_POS;
% % s = movmedian(POS.speed,100,'omitnan');
% % s = conv_filter(s,hanning(1500));
% % figure
% % plot(POS.Time_uS/60e6,s)
% % axis tight
% 
% %
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plot_me(OUT,LFP_files,sig_name)
% figure
% subplot(1,2,1)
% barh(OUT.INFO.ratings_non_reref)
% axis tight
% set(gca,'YTick',1:length(OUT.INFO.ratings_non_reref))
% set(gca,'YTickLabel',LFP_files,'FontSize',7)
% title(['Sig to Noise ' sig_name])
% subplot(1,2,2)
% imagesc(OUT.INFO.scores_by_combo)
% set(gca,'YTick',1:length(OUT.INFO.ratings_non_reref))
% set(gca,'YTickLabel',LFP_files,'FontSize',5)
% set(gca,'XTick',1:length(OUT.INFO.ratings_non_reref))
% set(gca,'XTickLabel',LFP_files,'FontSize',5)
% set(gca,'XTickLabelRotation',90)
% axis xy
% 
% end
% 
% function plot_all_LFP(L,LFP,LFP_files,sig_name)
% figure
% for ii = 1:Cols(L)
%     subplot(ceil(Cols(L)/8),8,ii)
%     pwelch(L(:,ii),LFP.sFreq,LFP.sFreq/2,256,LFP.sFreq)
%     title(LFP_files{ii},'FontSize',7)
% end
% sgtitle(sig_name)
% end
% 
% function plot_all_LFP_spec(LFP,LFP_files,sig_name)
% figure
% for iF = 1:length(LFP_files)
%     fqs = 1:.5:100;
%     LFP = LK_Load_and_Clean_LFP(LFP_files{iF});
%     
%     subplot(ceil(length(LFP_files)/8),8,iF)
%     spectrogram(LFP.LFP,LFP.sFreq*10,LFP.sFreq*5,fqs, LFP.sFreq,'yaxis')
%     %     pwelch(L(:,ii),LFP.sFreq,LFP.sFreq/2,256,LFP.sFreq)
%     title(LFP_files{iF},'FontSize',7)
% end
% sgtitle(sig_name)
% end