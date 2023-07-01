%%% Making logicals to pulls out conditions
IX_LID = categorical(TBL_M1.Group) == '6OHDA_LID';
IX_Lesion_M1 = IX_LID & categorical(TBL_M1.Hemisphere) == 'R';
IX_UnLesion_M1 = IX_LID & categorical(TBL_M1.Hemisphere) == 'L';
IX_SK = categorical(TBL_M1.Drugs) == 'Saline&Ketamine';
IX_Ket1_Les_M1 = IX_SK & IX_Lesion_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 1; % -20 to -5 min pre Ket inj if LDO&Ket day then this coincides with peak 80
IX_Ket1_UnLes_M1 = IX_SK & IX_UnLesion_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 1;
IX_Ket2_Les_M1 = IX_SK & IX_Lesion_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 2; % 2-30 min post Ket inj
IX_Ket2_UnLes_M1 = IX_SK & IX_UnLesion_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 2;
IX_Ket3_Les_M1 = IX_SK & IX_Lesion_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 3; % 60-80 min post Ket inj
IX_Ket3_UnLes_M1 = IX_SK & IX_UnLesion_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 3;

IX_LDO1_Les_M1 = IX_Lesion_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 1; % -20 to -5 min pre LDOPA inj
IX_LDO1_UnLes_M1 = IX_UnLesion_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 1;
IX_LDO2_Les_M1 = IX_Lesion_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 2; % 2-30 min post LDO inj
IX_LDO2_UnLes_M1 = IX_UnLesion_M1 &categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 2;
IX_LDO3_Les_M1 = IX_Lesion_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 3; % 30 - 50 min post LDOPA inj
IX_LDO3_UnLes_M1 = IX_UnLesion_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 3;
IX_LDO4_Les_M1 = IX_Lesion_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 4; % 60 - 80 min post LDOPA inj
IX_LDO4_UnLes_M1 = IX_UnLesion_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 4;
IX_LDO4_Les_M1_ket = IX_LDO4_Les_M1 & categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine'; % 60 - 80 min post LDOPA inj
IX_LDO4_UnLes_M1_ket = IX_LDO4_UnLes_M1 & categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine';
IX_LDO4_Les_M1_sal = IX_LDO4_Les_M1 & categorical(TBL_M1.Drugs) == 'LDOPA&Saline'; % 60 - 80 min post LDOPA inj
IX_LDO4_UnLes_M1_sal = IX_LDO4_UnLes_M1 & categorical(TBL_M1.Drugs) == 'LDOPA&Saline';

IX_CM1 = categorical(TBL_M1.Group) == 'SHAM';
IX_Right_M1 = IX_CM1 & categorical(TBL_M1.Hemisphere) == 'R';
IX_Left_M1 = IX_CM1 & categorical(TBL_M1.Hemisphere) == 'L';

IX_LDO1_CM1_R = IX_Right_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 1; % -20 to -5 min pre Ket inj if LDO&Ket day then this coincides with peak 80
IX_LDO2_CM1_R = IX_Right_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 2; % 2-30 min post Ket inj
IX_LDO3_CM1_R = IX_Right_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 3; % 60-80 min post Ket inj
IX_LDO4_CM1_R = IX_Right_M1 & categorical(TBL_M1.Condition) == 'LDO' & TBL_M1.Interval == 4; % 60-80 min post Ket inj
IX_LDO4_CM1R_ket = IX_LDO4_CM1_R & categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine'; % 60 - 80 min post LDOPA inj
IX_LDO4_CM1R_sal = IX_LDO4_CM1_R & categorical(TBL_M1.Drugs) == 'LDOPA&Saline'; % 60 - 80 min post LDOPA inj

IX_Ket1_CM1_R = IX_SK & IX_Right_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 1; % -20 to -5 min pre Ket inj if LDO&Ket day then this coincides with peak 80
IX_Ket2_CM1_R = IX_SK & IX_Right_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 2; % 2-30 min post Ket inj
IX_Ket3_CM1_R = IX_SK & IX_Right_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 3; % 60-80 min post Ket inj

IX_Ket1_CM1_L = IX_SK & IX_Left_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 1; % -20 to -5 min pre Ket inj if LDO&Ket day then this coincides with peak 80
IX_Ket2_CM1_L = IX_SK & IX_Left_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 2; % 2-30 min post Ket inj
IX_Ket3_CM1_L = IX_SK & IX_Left_M1 & categorical(TBL_M1.Condition) == 'KET' & TBL_M1.Interval == 3; % 60-80 min post Ket inj


% IX_Ket2_Les_CM1 = IX_Ket2_CM1 & categorical(TBL_M1.Hemisphere) == 'R' ;

%% Radial historgram plots
% Phase locking to beta during baseline


rad_binsize = deg2rad(15);
rad_edges = -pi:rad_binsize:pi;

%% LID animals LDO&Ket and LDO&Sal sessions
TBL_LID_Base_Les = TBL_M1(IX_LDO1_Les_M1,:);

fq = [1 2 3 4 5];
clr = lines(5);

 for iN = 1:size(TBL_LID_Base_Les)
     if TBL_LID_Base_Les.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Base_Les.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Base_Les.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_LID_Base_Les.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Base_Les.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_Base_Les.Rat(iN),iN,TBL_LID_Base_Les.Ang_p(iN,fq),TBL_LID_Base_Les.Ang_z(iN,fq),TBL_LID_Base_Les.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 for ii = 1: length(fq)
     figure
     histogram(TBL_LID_Base_Les.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g Baseline period Lesioned M1 LID n = %d sig n = %h', fq(ii), length(TBL_LID_Base_Les.Ang_to_shuf_p(:,ii)), sum(TBL_LID_Base_Les.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
 end
 
TBL_LID_LDO1_Les = TBL_M1(IX_LDO2_Les_M1,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_LID_LDO1_Les.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g immediately post LDO period Lesioned M1 LID n = %d sig n = %h', fq(ii), length(TBL_LID_LDO1_Les.Ang_to_shuf_p(:,ii)), sum(TBL_LID_LDO1_Les.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end

for iN = 1:size(TBL_LID_LDO1_Les)
     if TBL_LID_LDO1_Les.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_LDO1_Les.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_LDO1_Les.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_LID_LDO1_Les.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_LDO1_Les.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_LDO1_Les.Rat(iN),iN,TBL_LID_LDO1_Les.Ang_p(iN,fq),TBL_LID_LDO1_Les.Ang_z(iN,fq),TBL_LID_LDO1_Les.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 
TBL_LID_LDO2_Les = TBL_M1(IX_LDO3_Les_M1,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_LID_LDO2_Les.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g peak LDO period Lesioned M1 LID n = %d sig n = %h', fq(ii), length(TBL_LID_LDO2_Les.Ang_to_shuf_p(:,ii)), sum(TBL_LID_LDO2_Les.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end
 
TBL_LID_LDOKet_Les = TBL_M1(IX_LDO4_Les_M1_ket,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_LID_LDOKet_Les.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g post ket period Lesioned M1 LID n = %d sig n = %h', fq(ii), length(TBL_LID_LDOKet_Les.Ang_to_shuf_p(:,ii)), sum(TBL_LID_LDOKet_Les.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end
 
TBL_LID_LDOSal_Les = TBL_M1(IX_LDO4_Les_M1_sal,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_LID_LDOSal_Les.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g post saline period Lesioned M1 LID n = %d sig n = %h', fq(ii), length(TBL_LID_LDOSal_Les.Ang_to_shuf_p(:,ii)), sum(TBL_LID_LDOSal_Les.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end
 
% Saline and ketamine sessions LID animals
TBL_LID_Base_Les_SK = TBL_M1(IX_Ket1_Les_M1,:);

fq = [1 2 3 4 5];
clr = lines(5);

 for iN = 1:size(TBL_LID_Base_Les_SK)
     if TBL_LID_Base_Les_SK.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Base_Les_SK.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Base_Les_SK.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_LID_Base_Les_SK.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Base_Les_SK.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_Base_Les_SK.Rat(iN),iN,TBL_LID_Base_Les_SK.Ang_p(iN,fq),TBL_LID_Base_Les_SK.Ang_z(iN,fq),TBL_LID_Base_Les_SK.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 for ii = 1: length(fq)
     figure
     histogram(TBL_LID_Base_Les_SK.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g Baseline period Lesioned M1 LID n = %d sig n = %h', fq(ii), length(TBL_LID_Base_Les_SK.Ang_to_shuf_p(:,ii)), sum(TBL_LID_Base_Les_SK.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
 end
 
TBL_LID_Ket1_Les_SK = TBL_M1(IX_Ket2_Les_M1,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_LID_Ket1_Les_SK.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g immediately post ketamine period Lesioned M1 LID Sal&Ket n = %d sig n = %h', fq(ii), length(TBL_LID_Ket1_Les_SK.Ang_to_shuf_p(:,ii)), sum(TBL_LID_Ket1_Les_SK.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end

for iN = 1:size(TBL_LID_Ket1_Les_SK)
     if TBL_LID_Ket1_Les_SK.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Ket1_Les_SK.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Ket1_Les_SK.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_LID_Ket1_Les_SK.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Ket1_Les_SK.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_Ket1_Les_SK.Rat(iN),iN,TBL_LID_Ket1_Les_SK.Ang_p(iN,fq),TBL_LID_Ket1_Les_SK.Ang_z(iN,fq),TBL_LID_Ket1_Les_SK.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 
TBL_LID_Ket2_Les_SK = TBL_M1(IX_Ket3_Les_M1,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_LID_Ket2_Les_SK.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g 60 min post ket period Lesioned M1 LID SAl&Ket n = %d sig n = %h', fq(ii), length(TBL_LID_Ket2_Les_SK.Ang_to_shuf_p(:,ii)), sum(TBL_LID_Ket2_Les_SK.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end

for iN = 1:size(TBL_LID_Ket2_Les_SK)
     if TBL_LID_Ket2_Les_SK.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Ket2_Les_SK.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Ket2_Les_SK.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_LID_Ket2_Les_SK.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Ket2_Les_SK.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_Ket2_Les_SK.Rat(iN),iN,TBL_LID_Ket2_Les_SK.Ang_p(iN,fq),TBL_LID_Ket2_Les_SK.Ang_z(iN,fq),TBL_LID_Ket2_Les_SK.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
end

figure
bar([sum(TBL_LID_Base_Les.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_LID_LDO1_Les.Ang_to_shuf_p(:,fq) < 0.05) ...
    sum(TBL_LID_LDO2_Les.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_LID_LDOKet_Les.Ang_to_shuf_p(:,fq) < 0.05)])
pubify_figure_axis
title('LID phaselocking theta Ldo&ket')

figure
bar([sum(TBL_LID_Base_Les_SK.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_LID_Ket1_Les_SK.Ang_to_shuf_p(:,fq) < 0.05) ...
    sum(TBL_LID_Ket2_Les_SK.Ang_to_shuf_p(:,fq) < 0.05)])
pubify_figure_axis
title('LID phaselocking theta sal&ket')

fq = 3;
figure
bar([sum(TBL_LID_Base_Les.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_LID_LDO1_Les.Ang_to_shuf_p(:,fq) < 0.05) ...
    sum(TBL_LID_LDO2_Les.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_LID_LDOKet_Les.Ang_to_shuf_p(:,fq) < 0.05)])
pubify_figure_axis
title('LID phaselocking beta Ldo&ket')

figure
bar([sum(TBL_LID_Base_Les_SK.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_LID_Ket1_Les_SK.Ang_to_shuf_p(:,fq) < 0.05) ...
    sum(TBL_LID_Ket2_Les_SK.Ang_to_shuf_p(:,fq) < 0.05)])
pubify_figure_axis
title('LID phaselocking beta sal&ket')
%% SHAM animals LDO&Ket and LDO&Sal sessions
TBL_CM1R_Base = TBL_M1(IX_LDO1_CM1_R,:);

fq = [1 2 3 4 5];
clr = lines(5);

 for iN = 1:size(TBL_CM1R_Base)
     if TBL_CM1R_Base.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Base.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Base.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_CM1R_Base.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Base.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_CM1R_Base.Rat(iN),iN,TBL_CM1R_Base.Ang_p(iN,fq),TBL_CM1R_Base.Ang_z(iN,fq),TBL_CM1R_Base.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 for ii = 1: length(fq)
     figure
     histogram(TBL_CM1R_Base.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g Baseline period SHAM M1 right hem n = %d sig n = %h', fq(ii), length(TBL_CM1R_Base.Ang_to_shuf_p(:,ii)), sum(TBL_CM1R_Base.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
 end
 
TBL_CM1R_LDO1 = TBL_M1(IX_LDO2_CM1_R,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_CM1R_LDO1.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g immediately post LDO period SHAM M1 right hem n = %d sig n = %h', fq(ii), length(TBL_CM1R_LDO1.Ang_to_shuf_p(:,ii)), sum(TBL_CM1R_LDO1.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end

for iN = 1:size(TBL_CM1R_LDO1)
     if TBL_CM1R_LDO1.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_LDO1.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_LDO1.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_CM1R_LDO1.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_LDO1.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_CM1R_LDO1.Rat(iN),iN,TBL_CM1R_LDO1.Ang_p(iN,fq),TBL_CM1R_LDO1.Ang_z(iN,fq),TBL_CM1R_LDO1.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 
TBL_CM1R_LDO2 = TBL_M1(IX_LDO3_CM1_R,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_CM1R_LDO2.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g peak LDO period SHAM M1 right hem n = %d sig n = %h', fq(ii), length(TBL_CM1R_LDO2.Ang_to_shuf_p(:,ii)), sum(TBL_CM1R_LDO2.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end

for iN = 1:size(TBL_CM1R_LDO2)
     if TBL_CM1R_LDO2.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_LDO2.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_LDO2.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_CM1R_LDO2.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_LDO2.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_CM1R_LDO2.Rat(iN),iN,TBL_CM1R_LDO2.Ang_p(iN,fq),TBL_CM1R_LDO2.Ang_z(iN,fq),TBL_CM1R_LDO2.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
TBL_CM1R_LDOKet = TBL_M1(IX_LDO4_CM1R_ket,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_CM1R_LDOKet.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g post ket period SHAM M1 right hem n = %d sig n = %h', fq(ii), length(TBL_CM1R_LDOKet.Ang_to_shuf_p(:,ii)), sum(TBL_CM1R_LDOKet.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end
 
TBL_CM1R_LDOSal = TBL_M1(IX_LDO4_CM1R_sal,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_CM1R_LDOSal.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g post saline period SHAM M1 right hem n = %d sig n = %h', fq(ii), length(TBL_CM1R_LDOSal.Ang_to_shuf_p(:,ii)), sum(TBL_CM1R_LDOSal.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end
 
% SHAM saline and ketamine sessions 

TBL_CM1R_Base_SK = TBL_M1(IX_Ket1_CM1_R,:);

fq = [1 2 3 4 5];
clr = lines(5);

 for iN = 1:size(TBL_CM1R_Base_SK)
     if TBL_CM1R_Base_SK.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Base_SK.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Base_SK.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_CM1R_Base_SK.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Base_SK.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_CM1R_Base_SK.Rat(iN),iN,TBL_CM1R_Base_SK.Ang_p(iN,fq),TBL_CM1R_Base_SK.Ang_z(iN,fq),TBL_CM1R_Base_SK.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 for ii = 1: length(fq)
     figure
     histogram(TBL_CM1R_Base_SK.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g Baseline period Lesioned M1 LID n = %d sig n = %h', fq(ii), length(TBL_CM1R_Base_SK.Ang_to_shuf_p(:,ii)), sum(TBL_CM1R_Base_SK.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
 end
 
TBL_CM1R_Ket1_SK = TBL_M1(IX_Ket2_CM1_R,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_CM1R_Ket1_SK.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g immediately post ketamine period Lesioned M1 LID Sal&Ket n = %d sig n = %h', fq(ii), length(TBL_CM1R_Ket1_SK.Ang_to_shuf_p(:,ii)), sum(TBL_CM1R_Ket1_SK.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end

for iN = 1:size(TBL_CM1R_Ket1_SK)
     if TBL_CM1R_Ket1_SK.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Ket1_SK.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Ket1_SK.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_CM1R_Ket1_SK.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Ket1_SK.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_CM1R_Ket1_SK.Rat(iN),iN,TBL_CM1R_Ket1_SK.Ang_p(iN,fq),TBL_CM1R_Ket1_SK.Ang_z(iN,fq),TBL_CM1R_Ket1_SK.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 
TBL_CM1R_Ket2_SK = TBL_M1(IX_Ket3_CM1_R,:);

for ii = 1: length(fq)
     figure
     histogram(TBL_CM1R_Ket2_SK.Ang_to_shuf_p(:,ii),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', clr(ii,:) ) %[0.92 0.04 0.04]
     plot_ref_line(0.05,'line_width',2,'style',':')
     title(sprintf('Frequency %g 60 min post ket period Lesioned M1 LID SAl&Ket n = %d sig n = %h', fq(ii), length(TBL_CM1R_Ket2_SK.Ang_to_shuf_p(:,ii)), sum(TBL_CM1R_Ket2_SK.Ang_to_shuf_p(:,ii) < 0.05)),'FontSize',10)
end

for iN = 1:size(TBL_CM1R_Ket2_SK)
     if TBL_CM1R_Ket2_SK.Ang_to_shuf_p(iN,fq) < 0.05
%      if TBL_LID_Beta_Les.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Ket2_SK.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Ket2_SK.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_CM1R_Ket2_SK.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1R_Ket2_SK.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_CM1R_Ket2_SK.Rat(iN),iN,TBL_CM1R_Ket2_SK.Ang_p(iN,fq),TBL_CM1R_Ket2_SK.Ang_z(iN,fq),TBL_CM1R_Ket2_SK.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
end
 
figure
bar([sum(TBL_CM1R_Base.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_CM1R_LDO1.Ang_to_shuf_p(:,fq) < 0.05) ...
    sum(TBL_CM1R_LDO2.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_CM1R_LDOKet.Ang_to_shuf_p(:,fq) < 0.05)])
pubify_figure_axis
title('SHAM phaselocking theta Ldo&ket')

figure
bar([sum(TBL_CM1R_Base_SK.Ang_to_shuf_p(:,fq) < 0.05) sum(TBL_CM1R_Ket1_SK.Ang_to_shuf_p(:,fq) < 0.05) ...
    sum(TBL_CM1R_Ket2_SK.Ang_to_shuf_p(:,fq) < 0.05)])
pubify_figure_axis
title('SHAM phaselocking theta sal&ket')


%% 
TBL_LID_Beta_UnLes = TBL_M1(IX_LDO1_UnLes_M1,:);

fq = 3;


 for iN = 1:size(TBL_LID_Beta_UnLes)
     if TBL_LID_Beta_UnLes.Ang_p(iN,fq) < 0.05 && TBL_LID_Beta_UnLes.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Beta_UnLes.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Beta_UnLes.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_LID_Beta_UnLes.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_Beta_UnLes.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_Beta_UnLes.Rat(iN),iN,TBL_LID_Beta_UnLes.Ang_p(iN,fq),TBL_LID_Beta_UnLes.Ang_z(iN,fq),TBL_LID_Beta_UnLes.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
 
 figure
 histogram(TBL_LID_Beta_UnLes.Ang_to_shuf_p(:,2),40,'DisplayStyle', 'stairs', 'LineWidth', 2.5,'EdgeColor', [0.92 0.04 0.04])
% phase locking to 80Hz during peak 80Hz time in the lesioned hemisphere of LID
% rats



TBL_LID_80_Les = TBL_M1(IX_LDO3_Les_M1,:);

 for iN = 1:size(TBL_LID_80_Les)
        if TBL_LID_80_Les.Ang_p(iN,4) < 0.05
            figure(iN)
            clf
            bar(TBL_LID_80_Les.Ang_z(iN,4))
            set(gca,'XTickLabel',TBL_LID_80_Les.fq_ctrs(iN,4))
            ylabel('Rayleigh z')
            yyaxis right
            stem(TBL_LID_80_Les.Ang_p(iN,4))
            ylabel('p value')
            xlabel('Hz')
            title(sprintf('NEURON %g Rat %d',iN,TBL_LID_80_Les.Rat(iN)))
        end
 end

 for iN = 1:size(TBL_LID_80_Les)
      if TBL_LID_80_Les.Ang_p(iN,4) < 0.05 && TBL_LID_80_Les.Ang_to_shuf_p(iN,4) < 0.05
         figure
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_80_Les.sh_hist_rad_mn{iN,1}(4,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_80_Les.sh_hist_rad_mn{iN,1}(4,:)  + TBL_LID_80_Les.sh_hist_rad_95ci{iN,1}(4,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_80_Les.hist_rad{iN,1}(4,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_80_Les.Rat(iN),iN,TBL_LID_80_Les.Ang_p(iN,4),TBL_LID_80_Les.Ang_z(iN,4),TBL_LID_80_Les.Ang_to_shuf_p(iN,4)),'FontSize',10)
     end
 end
Sig_PL = find(TBL_LID_80_Les.Ang_p(:,4) < 0.05); % 43 sig out of 207
Les_80_IN = find(TBL_LID_80_Les.Neuron_type == 1 & TBL_LID_80_Les.Ang_p(:,4) < 0.05); % 3/43
Les_80_PY = find(TBL_LID_80_Les.Neuron_type == 2 & TBL_LID_80_Les.Ang_p(:,4) < 0.05); % 14/43
Les_80_PY100 = find(TBL_LID_80_Les.Neuron_type == 3 & TBL_LID_80_Les.Ang_p(:,4) < 0.05); % 26/43
Les_IN = find(TBL_LID_80_Les.Neuron_type == 1); % 14 21% of all IN are sig to 80
Les_PY = find(TBL_LID_80_Les.Neuron_type == 2); % 106 13% of all PY are sig to 80
Les_PY100 = find(TBL_LID_80_Les.Neuron_type == 3); % 87 29% of all PY100 are sig to 80
Sig_Pl_80_Les = TBL_LID_80_Les(TBL_LID_80_Les.Ang_p(:,4) < 0.05, :);

%% Rat 337 Session 7 deep dive on teh sig phase locked neuron
test = TBL_M1.Rat == 337 & TBL_M1.Session == 7 & TBL_M1.Depth_uM == 1.093750000000000e+03;
Rat337Ses7 = TBL_M1(test,:);

for iN = 1:size(Rat337Ses7)
    figure
    %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.sh_hist_rad_mn{iN,1}(4,:) ,'DisplayStyle','stairs','EdgeColor','k')
    hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.sh_hist_rad_mn{iN,1}(4,:)  + Rat337Ses7.sh_hist_rad_95ci{iN,1}(4,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
    %hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.hist_rad{iN,1}(4,:) ,'FaceAlpha',.5)
    title(sprintf('Rat %g Condition %s\n\r Interval %g p=%1.3f z=%2.1f pshuff=%1.3f to 80Hz',Rat337Ses7.Rat(iN),categorical(Rat337Ses7.Condition(iN)),Rat337Ses7.Interval(iN),Rat337Ses7.Ang_p(iN,4),Rat337Ses7.Ang_z(iN,4),Rat337Ses7.Ang_to_shuf_p(iN,4)),'FontSize',12)
end

for iN = 1:size(Rat337Ses7)
    figure
    %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.sh_hist_rad_mn{iN,1}(3,:) ,'DisplayStyle','stairs','EdgeColor','k')
    hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.sh_hist_rad_mn{iN,1}(3,:)  + Rat337Ses7.sh_hist_rad_95ci{iN,1}(3,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
    %hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.hist_rad{iN,1}(3,:) ,'FaceAlpha',.5)
    title(sprintf('Rat %g Condition %s\n\r Interval %g p=%1.3f z=%2.1f pshuff=%1.3f to 50Hz',Rat337Ses7.Rat(iN),categorical(Rat337Ses7.Condition(iN)),Rat337Ses7.Interval(iN),Rat337Ses7.Ang_p(iN,3),Rat337Ses7.Ang_z(iN,3),Rat337Ses7.Ang_to_shuf_p(iN,3)),'FontSize',12)
end 

for iN = 1:size(Rat337Ses7)
    figure
    %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
    hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.sh_hist_rad_mn{iN,1}(fq,:)  + Rat337Ses7.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
    %hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat337Ses7.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
    title(sprintf('Rat %g Condition %s\n\r Interval %g p=%1.3f z=%2.1f pshuff=%1.3f to beta Hz',Rat337Ses7.Rat(iN),categorical(Rat337Ses7.Condition(iN)),Rat337Ses7.Interval(iN),Rat337Ses7.Ang_p(iN,1),Rat337Ses7.Ang_z(iN,1),Rat337Ses7.Ang_to_shuf_p(iN,1)),'FontSize',12)
end  

for iN = 1:size(Rat320Ses20)
    figure
    %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat320Ses20.sh_hist_rad_mn{iN,1}(4,:) ,'DisplayStyle','stairs','EdgeColor','k')
    hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat320Ses20.sh_hist_rad_mn{iN,1}(4,:)  + Rat320Ses20.sh_hist_rad_95ci{iN,1}(4,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
    %hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat320Ses20.hist_rad{iN,1}(4,:) ,'FaceAlpha',.5)
    title(sprintf('Rat %g Condition %s\n\r Interval %g p=%1.3f z=%2.1f pshuff=%1.3f to 80Hz',Rat320Ses20.Rat(iN),categorical(Rat320Ses20.Condition(iN)),Rat320Ses20.Interval(iN),Rat320Ses20.Ang_p(iN,4),Rat320Ses20.Ang_z(iN,4),Rat320Ses20.Ang_to_shuf_p(iN,4)),'FontSize',12)
end

for iN = 1:size(Rat320Ses20)
    figure
    %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat320Ses20.sh_hist_rad_mn{iN,1}(3,:) ,'DisplayStyle','stairs','EdgeColor','k')
    hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat320Ses20.sh_hist_rad_mn{iN,1}(3,:)  + Rat320Ses20.sh_hist_rad_95ci{iN,1}(3,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
    %hold on
    polarhistogram('BinEdges',rad_edges,'BinCounts',Rat320Ses20.hist_rad{iN,1}(3,:) ,'FaceAlpha',.5)
    title(sprintf('Rat %g Condition %s\n\r Interval %g p=%1.3f z=%2.1f pshuff=%1.3f to 50Hz',Rat320Ses20.Rat(iN),categorical(Rat320Ses20.Condition(iN)),Rat320Ses20.Interval(iN),Rat320Ses20.Ang_p(iN,3),Rat320Ses20.Ang_z(iN,3),Rat320Ses20.Ang_to_shuf_p(iN,3)),'FontSize',12)
end 

%%
for test = 1:size(Sig_Pl_80_Les)
    figure
    plot(Sig_Pl_80_Les.AC{test,1})
end

for iN = 1:size(Sig_Pl_80_Les)
    Hist_rad_80(iN,:) = Sig_Pl_80_Les.hist_rad{iN,1}(4,:);
    Sh_Hist_rad_80(iN,:) = Sig_Pl_80_Les.sh_hist_rad_mn{iN,1}(4,:);
end

for iN = 1:size(TBL_LID_80_Les)
    All_Hist_rad_80(iN,:) = TBL_LID_80_Les.hist_rad{iN,1}(4,:);
    All_Sh_Hist_rad_80(iN,:) = TBL_LID_80_Les.sh_hist_rad_mn{iN,1}(4,:);
end

figure
histogram('BinEdges',rad_edges,'BinCounts',nanmean(All_Hist_rad_80))
hold on
histogram('BinEdges',rad_edges,'BinCounts',nanmean(All_Sh_Hist_rad_80))
title('Sig phase locking 80Hz LID lesioned hemisphere mean hist rad with shuffeled data All neurons regardless of sig')

% polar histograms of the mean hist rad edges
figure
polarhistogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_80),'DisplayStyle','stairs','EdgeColor','k')
hold on
polarhistogram('BinEdges',rad_edges,'BinCounts',mean(Sh_Hist_rad_80), 'FaceAlpha',.5 )
title('Sig phase locking 80Hz LID lesioned hemisphere mean hist rad with shuffeled data polar plots')




% effect size calculation; cohens D
Hist_rad_80_D = Cohens_d(Hist_rad_80);
Ang_z_80_D = Cohens_d(Sig_Pl_80_Les.Ang_z(:,4));

for iN = 1:length(Les_80_IN)
    in = Les_80_IN(iN);
    Hist_rad_80_IN(iN,:) = TBL_LID_80_Les.hist_rad{in,1}(4,:);
    
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_80_IN))

for iN = 1:length(Les_80_PY)
    in = Les_80_PY(iN);
    Hist_rad_80_PY(iN,:) = TBL_LID_80_Les.hist_rad{in,1}(4,:);
    
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_80_PY))


for iN = 1:length(Les_80_PY100)
    in = Les_80_PY100(iN);
    Hist_rad_80_PY100(iN,:) = TBL_LID_80_Les.hist_rad{in,1}(4,:);
    
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_80_PY100))


for iu = 1:length(TS)
    figure;
    [PhaseLocking_80,ix,x_sec] = PETH_EEG_simple([LFP.t_uS(ix1:ix2), L(:,4)],TSr{1,iu},LFP.sFreq/50,LFP.sFreq/50,LFP.sFreq,true);
    colorbar
    title(sprintf ('Neuron %0.1f Rat 342 Ses 6 Peak 80Hz period', iu))
end

%Just checking the UnLes hemisphere
TBL_LID_80_UnLes = TBL_M1(IX_LDO3_UnLes_M1,:);

Sig_Pl_80_UnLes = TBL_LID_80_UnLes(TBL_LID_80_UnLes.Ang_p(:,4) < 0.05, :);
for iN = 1:size(Sig_Pl_80_UnLes)
    Hist_rad_80_UnLes(iN,:) = Sig_Pl_80_UnLes.hist_rad{iN,1}(4,:);
    Sh_Hist_rad_80_UnLes(iN,:) = Sig_Pl_80_UnLes.sh_hist_rad_mn{iN,1}(4,:);
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_80_UnLes))
hold on
histogram('BinEdges',rad_edges,'BinCounts',mean(Sh_Hist_rad_80_UnLes))
title('Sig phase locking 80Hz LID Unlesioned hemisphere mean hist rad with shuffeled data')

for iN = 1:size(TBL_LID_80_UnLes)
      if TBL_LID_80_UnLes.Ang_p(iN,4) < 0.05 && TBL_LID_80_UnLes.Ang_to_shuf_p(iN,4) < 0.05
         figure
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_80_UnLes.sh_hist_rad_mn{iN,1}(4,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_80_UnLes.sh_hist_rad_mn{iN,1}(4,:)  + TBL_LID_80_UnLes.sh_hist_rad_95ci{iN,1}(4,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_80_UnLes.hist_rad{iN,1}(4,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_80_UnLes.Rat(iN),iN,TBL_LID_80_UnLes.Ang_p(iN,4),TBL_LID_80_UnLes.Ang_z(iN,4),TBL_LID_80_UnLes.Ang_to_shuf_p(iN,4)),'FontSize',10)
     end
end

 
%% phase locking for 80 during the first 20 min post LDOPA inj just to see whats there
TBL_LID_No80_Les = TBL_M1(IX_LDO2_Les_M1,:);

 for iN = 1:size(TBL_LID_No80_Les)
      if TBL_LID_No80_Les.Ang_p(iN,4) < 0.05 && TBL_LID_No80_Les.Ang_to_shuf_p(iN,4) < 0.05
         figure
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_No80_Les.sh_hist_rad_mn{iN,1}(4,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_No80_Les.sh_hist_rad_mn{iN,1}(4,:)  + TBL_LID_No80_Les.sh_hist_rad_95ci{iN,1}(4,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_No80_Les.hist_rad{iN,1}(4,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_No80_Les.Rat(iN),iN,TBL_LID_No80_Les.Ang_p(iN,4),TBL_LID_No80_Les.Ang_z(iN,4),TBL_LID_No80_Les.Ang_to_shuf_p(iN,4)),'FontSize',10)
     end
 end

 
TBL_LID_BeforeLDO_Les = TBL_M1(IX_LDO1_Les_M1,:);

 for iN = 1:size(TBL_LID_BeforeLDO_Les)
      if TBL_LID_BeforeLDO_Les.Ang_p(iN,4) < 0.05 && TBL_LID_BeforeLDO_Les.Ang_to_shuf_p(iN,4) < 0.05
         figure
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_BeforeLDO_Les.sh_hist_rad_mn{iN,1}(4,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_BeforeLDO_Les.sh_hist_rad_mn{iN,1}(4,:)  + TBL_LID_BeforeLDO_Les.sh_hist_rad_95ci{iN,1}(4,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_BeforeLDO_Les.hist_rad{iN,1}(4,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_BeforeLDO_Les.Rat(iN),iN,TBL_LID_BeforeLDO_Les.Ang_p(iN,4),TBL_LID_BeforeLDO_Les.Ang_z(iN,4),TBL_LID_BeforeLDO_Les.Ang_to_shuf_p(iN,4)),'FontSize',10)
     end
 end
 

%% phase locking to 50Hz during peak Ketamine time in the lesioned hemisphere of LID
% rats
TBL_LID_50_Les = TBL_M1(IX_Ket2_Les_M1,:);

 for iN = 1:size(TBL_LID_50_Les)
        if TBL_LID_50_Les.Ang_p(iN,3) < 0.05
            figure(iN)
            clf
            bar(TBL_LID_50_Les.Ang_z(iN,3))
            set(gca,'XTickLabel',TBL_LID_50_Les.fq_ctrs(iN,3))
            ylabel('Rayleigh z')
            yyaxis right
            stem(TBL_LID_50_Les.Ang_p(iN,3))
            ylabel('p value')
            xlabel('Hz')
            title(sprintf('NEURON %g Rat %d',iN,TBL_LID_50_Les.Rat(iN)))
        end
 end

 for iN = 1:size(TBL_LID_50_Les)
      if TBL_LID_50_Les.Ang_p(iN,3) < 0.05 && TBL_LID_50_Les.Ang_to_shuf_p(iN,3) > 0.05
         figure
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_Les.sh_hist_rad_mn{iN,1}(3,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_Les.sh_hist_rad_mn{iN,1}(3,:)  + TBL_LID_50_Les.sh_hist_rad_95ci{iN,1}(3,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_Les.hist_rad{iN,1}(3,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_50_Les.Rat(iN),iN,TBL_LID_50_Les.Ang_p(iN,3),TBL_LID_50_Les.Ang_z(iN,3),TBL_LID_50_Les.Ang_to_shuf_p(iN,3)),'FontSize',10)
     end
 end

 for iN = 1:size(TBL_LID_50_Les)
      if TBL_LID_50_Les.Ang_p(iN,1) < 0.05 && TBL_LID_50_Les.Ang_to_shuf_p(iN,1) > 0.05
         figure
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_Les.sh_hist_rad_mn{iN,1}(1,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_Les.sh_hist_rad_mn{iN,1}(1,:)  + TBL_LID_50_Les.sh_hist_rad_95ci{iN,1}(1,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_Les.hist_rad{iN,1}(1,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_50_Les.Rat(iN),iN,TBL_LID_50_Les.Ang_p(iN,1),TBL_LID_50_Les.Ang_z(iN,1),TBL_LID_50_Les.Ang_to_shuf_p(iN,1)),'FontSize',10)
     end
 end
 
Sig50_PL = find(TBL_LID_50_Les.Ang_p(:,3) < 0.05); % 65 sig out of 236
Sal_Ket_PL = find(categorical(TBL_LID_50_Les.Drugs) == 'Saline&Ketamine' & TBL_LID_50_Les.Ang_p(:,3) < 0.05); % 37/65
LDO_Ket_PL = find(categorical(TBL_LID_50_Les.Drugs) == 'LDOPA&Ketamine' & TBL_LID_50_Les.Ang_p(:,3) < 0.05); % 28/65

Les_50_IN = find(TBL_LID_50_Les.Neuron_type == 1 & TBL_LID_50_Les.Ang_p(:,3) < 0.05); % 4/65
Les_50_PY = find(TBL_LID_50_Les.Neuron_type == 2 & TBL_LID_50_Les.Ang_p(:,3) < 0.05); % 22/65
Les_50_PY100 = find(TBL_LID_50_Les.Neuron_type == 3 & TBL_LID_50_Les.Ang_p(:,3) < 0.05); % 39/65

Les_50_IN_SK = find(TBL_LID_50_Les.Neuron_type == 1 & TBL_LID_50_Les.Ang_p(:,3) < 0.05 & categorical(TBL_LID_50_Les.Drugs) == 'Saline&Ketamine'); % 4/37
Les_50_PY_SK = find(TBL_LID_50_Les.Neuron_type == 2 & TBL_LID_50_Les.Ang_p(:,3) < 0.05 & categorical(TBL_LID_50_Les.Drugs) == 'Saline&Ketamine'); % 12/37
Les_50_PY100_SK = find(TBL_LID_50_Les.Neuron_type == 3 & TBL_LID_50_Les.Ang_p(:,3) < 0.05 & categorical(TBL_LID_50_Les.Drugs) == 'Saline&Ketamine'); % 21/37

Les_50_IN_LK = find(TBL_LID_50_Les.Neuron_type == 1 & TBL_LID_50_Les.Ang_p(:,3) < 0.05 & categorical(TBL_LID_50_Les.Drugs) == 'LDOPA&Ketamine'); % 0/28
Les_50_PY_LK = find(TBL_LID_50_Les.Neuron_type == 2 & TBL_LID_50_Les.Ang_p(:,3) < 0.05 & categorical(TBL_LID_50_Les.Drugs) == 'LDOPA&Ketamine'); % 10/28
Les_50_PY100_LK = find(TBL_LID_50_Les.Neuron_type == 3 & TBL_LID_50_Les.Ang_p(:,3) < 0.05 & categorical(TBL_LID_50_Les.Drugs) == 'LDOPA&Ketamine'); % 18/28

Les50_IN = find(TBL_LID_50_Les.Neuron_type == 1); %16 25% 
Les50_PY = find(TBL_LID_50_Les.Neuron_type == 2); %112 19% 
Les50_PY100 = find(TBL_LID_50_Les.Neuron_type == 3); %104 37.5%

Sig_Pl_50_Les = TBL_LID_50_Les(TBL_LID_50_Les.Ang_p(:,3) < 0.05, :);
for iN = 1:size(Sig_Pl_50_Les)
    Hist_rad_50(iN,:) = Sig_Pl_50_Les.hist_rad{iN,1}(3,:);
    Sh_Hist_rad_50(iN,:) = Sig_Pl_50_Les.sh_hist_rad_mn{iN,1}(3,:);
    Ang_z_50(iN,:) = Sig_Pl_50_Les.Ang_z(iN,3);
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_50))
hold on
histogram('BinEdges',rad_edges,'BinCounts',mean(Sh_Hist_rad_50))
title('Sig phase locking 50Hz LID lesioned hemisphere mean hist rad with shuffeled data')

figure
polarhistogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_50),'DisplayStyle','stairs','EdgeColor','k')
hold on
polarhistogram('BinEdges',rad_edges,'BinCounts',mean(Sh_Hist_rad_50), 'FaceAlpha',.5 )
title('Sig phase locking 50Hz LID lesioned hemisphere mean hist rad with shuffeled data polar plots')


% Effect size calculation
Hist_rad_50_D = Cohens_d(Hist_rad_50);
Ang_z_50_D = Cohens_d(Ang_z_50);

for iN = 1:length(Les_50_IN)
    in = Les_50_IN(iN);
    Hist_rad_50_IN(iN,:) = TBL_LID_50_Les.hist_rad{in,1}(3,:);
    Ang_z_50_IN(iN,:) = TBL_LID_50_Les.Ang_z(iN,3);
    
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_50_IN))

Ang_z_50_D_IN = Cohens_d(Ang_z_50_IN);

for iN = 1:length(Les_50_PY)
    in = Les_50_PY(iN);
    Hist_rad_50_PY(iN,:) = TBL_LID_50_Les.hist_rad{in,1}(3,:);
    Sh_Hist_rad_50_PY(iN,:) = TBL_LID_50_Les.sh_hist_rad_mn{in,1}(3,:);
    Ang_z_50_PY(iN,:) = TBL_LID_50_Les.Ang_z(iN,3);
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_50_PY))
hold on
histogram('BinEdges',rad_edges,'BinCounts',mean(Sh_Hist_rad_50_PY))
title('Sig phase locking principal cells only to 50Hz LID lesioned hemisphere mean hist rad')

Ang_z_50_D_PY = Cohens_d(Ang_z_50_PY);

for iN = 1:length(Les_50_PY100)
    in = Les_50_PY100(iN);
    Hist_rad_50_PY100(iN,:) = TBL_LID_50_Les.hist_rad{in,1}(3,:);
    Ang_z_50_PY100(iN,:) = TBL_LID_50_Les.Ang_z(iN,3);
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_50_PY100))

Ang_z_50_D_PY100 = Cohens_d(Ang_z_50_PY100);

%% Phase locking to 50 Hz in unlesioned hemisphere 
TBL_LID_50_UnLes = TBL_M1(IX_Ket2_UnLes_M1,:);

 fq = 3;
 
 for iN = 1:size(TBL_LID_50_UnLes)
     if TBL_LID_50_UnLes.Ang_p(iN,fq) < 0.05 && TBL_LID_50_UnLes.Ang_to_shuf_p(iN,fq) < 0.05
         figure
         %          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_UnLes.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_UnLes.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_LID_50_UnLes.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_LID_50_UnLes.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_LID_50_UnLes.Rat(iN),iN,TBL_LID_50_UnLes.Ang_p(iN,fq),TBL_LID_50_UnLes.Ang_z(iN,fq),TBL_LID_50_UnLes.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end

%% Phase locking to 80Hz during peak LDO in sham lesioned hemisphere of
% shma rats

TBL_CM1_80_Les = TBL_M1(IX_LDO3_CM1_R,:);
fq = 1;

 for iN = 1:size(TBL_CM1_80_Les)
      if TBL_CM1_80_Les.Ang_p(iN,fq) < 0.05 && TBL_CM1_80_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure(iN)
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1_80_Les.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1_80_Les.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_CM1_80_Les.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1_80_Les.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_CM1_80_Les.Rat(iN),iN,TBL_CM1_80_Les.Ang_p(iN,fq),TBL_CM1_80_Les.Ang_z(iN,fq),TBL_CM1_80_Les.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end


%% Phase locking to 50Hz during peak ketamine in sham lesioned hemisphere of
% shma rats

TBL_CM1_50_Les = TBL_M1(IX_Ket2_CM1_R,:);

 for iN = 1:size(TBL_CM1_50_Les)
        if TBL_CM1_50_Les.Ang_p(iN,3) < 0.05
            figure(iN)
            clf
            bar(TBL_CM1_50_Les.Ang_z(iN,3))
            set(gca,'XTickLabel',TBL_CM1_50_Les.fq_ctrs(iN,3))
            ylabel('Rayleigh z')
            yyaxis right
            stem(TBL_CM1_50_Les.Ang_p(iN,3))
            ylabel('p value')
            xlabel('Hz')
            title(sprintf('NEURON %g Rat %d',iN,TBL_CM1_50_Les.Rat(iN)))
        end
 end

fq = 3;

 for iN = 1:size(TBL_CM1_50_Les)
      if TBL_CM1_50_Les.Ang_p(iN,fq) < 0.05 && TBL_CM1_50_Les.Ang_to_shuf_p(iN,fq) < 0.05
         figure
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1_50_Les.sh_hist_rad_mn{iN,1}(fq,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1_50_Les.sh_hist_rad_mn{iN,1}(fq,:)  + TBL_CM1_50_Les.sh_hist_rad_95ci{iN,1}(fq,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',TBL_CM1_50_Les.hist_rad{iN,1}(fq,:) ,'FaceAlpha',.5)
         title(sprintf('Rat %g Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',TBL_CM1_50_Les.Rat(iN),iN,TBL_CM1_50_Les.Ang_p(iN,fq),TBL_CM1_50_Les.Ang_z(iN,fq),TBL_CM1_50_Les.Ang_to_shuf_p(iN,fq)),'FontSize',10)
     end
 end
Sig50_PL_CM1 = find(TBL_CM1_50_Les.Ang_p(:,3) < 0.05); % 11 out of 87
Sal_Ket_CM1 = find(categorical(TBL_CM1_50_Les.Drugs) == 'Saline&Ketamine' & TBL_CM1_50_Les.Ang_p(:,3) < 0.05); % 4/9
LDO_Ket_CM1 = find(categorical(TBL_CM1_50_Les.Drugs) == 'LDOPA&Ketamine' & TBL_CM1_50_Les.Ang_p(:,3) < 0.05); % 7/11

CM1_50_IN = find(TBL_CM1_50_Les.Neuron_type == 1 & TBL_CM1_50_Les.Ang_p(:,3) < 0.05); % 2/11
CM1_50_PY = find(TBL_CM1_50_Les.Neuron_type == 2 & TBL_CM1_50_Les.Ang_p(:,3) < 0.05); % 3/11
CM1_50_PY100 = find(TBL_CM1_50_Les.Neuron_type == 3 & TBL_CM1_50_Les.Ang_p(:,3) < 0.05); % 6/11

for iN = 1:size(Sig50_PL_CM1)
    Hist_rad_CM_50(iN,:) = Sig50_PL_CM1.hist_rad{iN,1}(3,:);
    Sh_Hist_rad_CM_50(iN,:) = Sig50_PL_CM1.sh_hist_rad_mn{iN,1}(3,:);
    Ang_z_CM_50(iN,:) = Sig50_PL_CM1.Ang_z(iN,3);
end
figure
histogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_CM_50))
hold on
histogram('BinEdges',rad_edges,'BinCounts',mean(Sh_Hist_rad_CM_50))
title('Sig phase locking 50Hz LID lesioned hemisphere mean hist rad with shuffeled data')

figure
polarhistogram('BinEdges',rad_edges,'BinCounts',mean(Hist_rad_CM_50),'DisplayStyle','stairs','EdgeColor','k')
hold on
polarhistogram('BinEdges',rad_edges,'BinCounts',mean(Sh_Hist_rad_CM_50), 'FaceAlpha',.5 )
title('Sig phase locking 50Hz LID lesioned hemisphere mean hist rad with shuffeled data polar plots')


%% fill in empty Coh_z_center
Coh_z_center = nan(Rows(TBL),499);
Coh_p_center = nan(Rows(TBL),499);
for ii = 1:Rows(TBL)
    if ~isempty(TBL.Coh_circ_rtest_p_center{ii,1})
       Coh_z_center(ii,:) = TBL.Coh_circ_rtest_z_center{ii,1}' ;
       Coh_p_center(ii,:) = TBL.Coh_circ_rtest_p_center{ii,1}' ;
    end
end
%standardized
Coh_z_center_norm = standardize_range(Coh_z_center',[0 1])';
Coh_z_center_norm_sig = Coh_z_center_norm;
Coh_z_center_norm_sig(Coh_p_center > .05) = nan; % Delete the 1st and 488th and 499th columns
% IX = Coh_p_center > .05;
% Coh_z_center_norm_sig_only = Coh_z_center_norm(IX);
% Coh_z_center_norm_sig_ket1_les = Coh_z_center_norm(IX_Ket1_Les_M1,:);
% Coh_z_center_norm_sig_ket1_les = Coh_z_center_norm_sig_ket1_les(IX(IX_Ket1_Les_M1,:));
% Coh_z_center_norm_sig_ket2_les = Coh_z_center_norm(IX_Ket2_Les_M1,:);
% Coh_z_center_norm_sig_ket2_les = Coh_z_center_norm_sig_ket2_les(IX(IX_Ket2_Les_M1,:));
% Coh_z_center_norm_sig_ket3_les = Coh_z_center_norm(IX_Ket3_Les_M1,:);
% Coh_z_center_norm_sig_ket3_les = Coh_z_center_norm_sig_ket3_les(IX(IX_Ket3_Les_M1,:));
% Coh_z_center_norm_sig_LDO1_les = Coh_z_center_norm(IX_LDO1_Les_M1,:);
% Coh_z_center_norm_sig_LDO1_les = Coh_z_center_norm_sig_LDO1_les(IX(IX_LDO1_Les_M1,:));
% Correlate pre and post ketamine periods of z scores

%% Mean spike field coupling to plot confidence intervals on the imagesc plots below
% Lesioned hemisphere
% mn_Ket1_Les = nanmean(Coh_z_center_norm_sig(IX_Ket1_Les_M1,:));
% SEM
% mn_Ket2_Les = nanmean(Coh_z_center_norm_sig(IX_Ket2_Les_M1,:));
% mn_Ket3_Les = nanmean(Coh_z_center_norm_sig(IX_Ket3_Les_M1,:));
% mn_LDO1_Les = nanmean(Coh_z_center_norm_sig(IX_LDO1_Les_M1,:));
% mn_LDO3_Les = nanmean(Coh_z_center_norm_sig(IX_LDO3_Les_M1,:));
% mn = [];
% SEM = [];
% mn(1,:) = nanmean(Coh_z_center_norm_sig(IX_Ket1_Les_M1,:));
% SEM(1,:) = std(Coh_z_center_norm_sig(IX_Ket1_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm_sig(IX_Ket1_Les_M1,:)));
% mn(2,:) = nanmean(Coh_z_center_norm_sig(IX_Ket2_Les_M1,:));
% SEM(2,:) = std(Coh_z_center_norm_sig(IX_Ket2_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm_sig(IX_Ket2_Les_M1,:)));
% mn(3,:) = nanmean(Coh_z_center_norm_sig(IX_Ket3_Les_M1,:));
% SEM(3,:) = std(Coh_z_center_norm_sig(IX_Ket3_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm_sig(IX_Ket3_Les_M1,:)));
% mn(4,:) = nanmean(Coh_z_center_norm_sig(IX_LDO1_Les_M1,:));
% SEM(4,:) = std(Coh_z_center_norm_sig(IX_LDO1_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm_sig(IX_LDO1_Les_M1,:)));
% mn(5,:) = nanmean(Coh_z_center_norm_sig(IX_LDO3_Les_M1,:));
% SEM(5,:) = std(Coh_z_center_norm_sig(IX_LDO3_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm_sig(IX_LDO3_Les_M1,:)));

mn_Les = [];
SEM_Les = [];
mn_Les(1,:) = nanmean(Coh_z_center_norm(IX_Ket1_Les_M1,:));
SEM_Les(1,:) = std(Coh_z_center_norm(IX_Ket1_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_Ket1_Les_M1,:)));
mn_Les(2,:) = nanmean(Coh_z_center_norm(IX_Ket2_Les_M1,:));
SEM_Les(2,:) = std(Coh_z_center_norm(IX_Ket2_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_Ket2_Les_M1,:)));
mn_Les(3,:) = nanmean(Coh_z_center_norm(IX_Ket3_Les_M1,:));
SEM_Les(3,:) = std(Coh_z_center_norm(IX_Ket3_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_Ket3_Les_M1,:)));
mn_Les(4,:) = nanmean(Coh_z_center_norm(IX_LDO1_Les_M1,:));
SEM_Les(4,:) = std(Coh_z_center_norm(IX_LDO1_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_LDO1_Les_M1,:)));
mn_Les(5,:) = nanmean(Coh_z_center_norm(IX_LDO3_Les_M1,:));
SEM_Les(5,:) = std(Coh_z_center_norm(IX_LDO3_Les_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_LDO3_Les_M1,:)));

mn_UnLes = [];
SEM_UnLes = [];
mn_UnLes(1,:) = nanmean(Coh_z_center_norm(IX_Ket1_UnLes_M1,:));
SEM_UnLes(1,:) = std(Coh_z_center_norm(IX_Ket1_UnLes_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_Ket1_UnLes_M1,:)));
mn_UnLes(2,:) = nanmean(Coh_z_center_norm(IX_Ket2_UnLes_M1,:));
SEM_UnLes(2,:) = std(Coh_z_center_norm(IX_Ket2_UnLes_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_Ket2_UnLes_M1,:)));
mn_UnLes(3,:) = nanmean(Coh_z_center_norm(IX_Ket3_UnLes_M1,:));
SEM_UnLes(3,:) = std(Coh_z_center_norm(IX_Ket3_UnLes_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_Ket3_UnLes_M1,:)));
mn_UnLes(4,:) = nanmean(Coh_z_center_norm(IX_LDO1_UnLes_M1,:));
SEM_UnLes(4,:) = std(Coh_z_center_norm(IX_LDO1_UnLes_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_LDO1_UnLes_M1,:)));
mn_UnLes(5,:) = nanmean(Coh_z_center_norm(IX_LDO3_UnLes_M1,:));
SEM_UnLes(5,:) = std(Coh_z_center_norm(IX_LDO3_UnLes_M1,:),'omitnan')./sqrt(length(Coh_z_center_norm(IX_LDO3_UnLes_M1,:)));


psd_fqs = TBL.Coh_psd_fqs{1,1};
clrs = lines(5);
for ii = 1:length(mn(:,1))
    figure
    errorbar(psd_fqs(1:5:end),mn(ii,1:5:end),SEM(ii,1:5:end))
    title(ii)
    pubify_figure_axis
    ax=gca;
    ax.FontSize=35;
    set(gca,'ylim',[0 .95])
    set(gca,'xlim',[0 250])
end

% confidence intervals plotting
% Lesioned hemisphere
CohP_Les = [];
CohP_Les = Coh_z_center_norm(IX_LDO1_Les_M1,:);
CohP_Les = cat(3,CohP_Les,Coh_z_center_norm(IX_LDO3_Les_M1,:));
CohP_Les = cat(3,CohP_Les,Coh_z_center_norm(IX_Ket2_Les_M1,:));
CohP_Les = cat(3,CohP_Les,Coh_z_center_norm(IX_Ket3_Les_M1,:));
for ii = 1:length(CohP_Les(1,1,:))
    figure
    plot_confidence_intervals(TBL.Coh_psd_fqs{1,1},CohP_Les(:,:,ii),[],clrs(1,:))
    title(sprintf('Lesion %0.0f',ii))
    pubify_figure_axis
    ax=gca;
    ax.FontSize=35;
    set(gca,'ylim',[0 .95])
    xticks([0 100 200 300 400 500])
    xticklabels({'0','50','100','150','200','250'})
end

% Un-lesioned hemisphere
CohP = [];
CohP = Coh_z_center_norm(IX_LDO1_UnLes_M1,:);
CohP = cat(3,CohP,Coh_z_center_norm(IX_LDO3_UnLes_M1,:));
CohP = cat(3,CohP,Coh_z_center_norm(IX_Ket2_UnLes_M1,:));
CohP = cat(3,CohP,Coh_z_center_norm(IX_Ket3_UnLes_M1,:));
for ii = 1:length(CohP(1,1,:))
    figure
    plot_confidence_intervals(TBL.Coh_psd_fqs{1,1},CohP(:,:,ii),[],clrs(1,:))
    title(sprintf('UnLesion %0.0f',ii))
    pubify_figure_axis
    ax=gca;
    ax.FontSize=35;
    set(gca,'ylim',[0 .8])
    xticks([0 100 200 300 400 500])
    xticklabels({'0','50','100','150','200','250'})
end

% ANOVA
gamma_80 = [74 94];
IX_80 = psd_fqs > gamma_80(1) & psd_fqs < gamma_80(2);
gamma_60 = [60 90];
IX_60 = psd_fqs > gamma_60(1) & psd_fqs < gamma_60(2);
gamma_50 = [45 75];
IX_50 = psd_fqs > gamma_50(1) & psd_fqs < gamma_50(2);
% Lesioned hemisphere
mn_ANV_L(:,1) = nanmean(Coh_z_center_norm(IX_LDO1_Les_M1,IX_80)')';
mn_ANV_L(:,2) = nanmean(Coh_z_center_norm(IX_LDO3_Les_M1,IX_80)')';
mn_ANV_L(:,3) = nanmean(Coh_z_center_norm(IX_Ket2_Les_M1,IX_80)')';
mn_ANV_L(:,4) = nanmean(Coh_z_center_norm(IX_Ket3_Les_M1,IX_80)')';
[p_L,tbl,stats_l] = kruskalwallis(mn_ANV_L);
% [p_L,tbl,stats_l] = anova1(mn_ANV_L);
figure
pairwise_L = multcompare(stats_l);
%Unlesioned
mn_ANV(:,1) = nanmean(Coh_z_center_norm(IX_LDO1_UnLes_M1,IX_60)')';
mn_ANV(:,2) = nanmean(Coh_z_center_norm(IX_LDO3_UnLes_M1,IX_60)')';
mn_ANV(:,3) = nanmean(Coh_z_center_norm(IX_Ket2_UnLes_M1,IX_60)')';
mn_ANV(:,4) = nanmean(Coh_z_center_norm(IX_Ket3_UnLes_M1,IX_60)')';
[p_UL, tbl, stats] = kruskalwallis(mn_ANV);
% [p_UL,tbl,stats] = anova1(mn_ANV);
figure
pairwise = multcompare(stats);

% Rank sum 
p_L_rank = signrank(mn_ANV_L(:,2)-mn_ANV_L(:,3));

% Lesioned hemisphere
mn_ANV_L(:,1) = nanmean(Coh_z_center_norm(IX_LDO1_Les_M1,IX_60)')';
mn_ANV_L(:,2) = nanmean(Coh_z_center_norm(IX_LDO3_Les_M1,IX_60)')';
mn_ANV_L(:,3) = nanmean(Coh_z_center_norm(IX_Ket2_Les_M1,IX_60)')';
mn_ANV_L(:,4) = nanmean(Coh_z_center_norm(IX_Ket3_Les_M1,IX_60)')';
[p_L,tbl,stats_l] = kruskalwallis(mn_ANV_L);
pairwise_L = multcompare(stats_l);
%Unlesioned
mn_ANV(:,1) = nanmean(Coh_z_center_norm(IX_LDO1_UnLes_M1,IX_60)')';
mn_ANV(:,2) = nanmean(Coh_z_center_norm(IX_LDO3_UnLes_M1,IX_60)')';
mn_ANV(:,3) = nanmean(Coh_z_center_norm(IX_Ket2_UnLes_M1,IX_60)')';
mn_ANV(:,4) = nanmean(Coh_z_center_norm(IX_Ket3_UnLes_M1,IX_60)')';
[p_UL, tbl, stats] = kruskalwallis(mn_ANV);
pairwise = multcompare(stats);


% Control 
CohP_CM1 = [];
CohP_CM1 = Coh_z_center_norm(IX_Ket1_CM1,:);
CohP_CM1 = cat(3,CohP_CM1,Coh_z_center_norm(IX_Ket2_CM1,:));
CohP_CM1 = cat(3,CohP_CM1,Coh_z_center_norm(IX_Ket3_CM1,:));
for ii = 1:length(CohP_CM1(1,1,:))
    figure
    plot_confidence_intervals(TBL.Coh_psd_fqs{1,1},CohP_CM1(:,:,ii),[],clrs(1,:))
    title(sprintf('Control %0.0f',ii))
    pubify_figure_axis
    ax=gca;
    ax.FontSize=35;
    set(gca,'ylim',[.1 .5])
%     xticks([0 100 200 300 400 500])
%     xticklabels({'0','50','100','150','200','250'})
end

% ANOVA
mn_ANV_CM1(:,1) = nanmean(Coh_z_center_norm(IX_Ket1_CM1,IX_50)')';
mn_ANV_CM1(:,2) = nanmean(Coh_z_center_norm(IX_Ket2_CM1,IX_50)')';
mn_ANV_CM1(:,3) = nanmean(Coh_z_center_norm(IX_Ket3_CM1,IX_50)')';

[p_CM1, tbl, stats_CM1] = kruskalwallis(mn_ANV_CM1);
figure
pairwise_CM1 = multcompare(stats_CM1);

%% Correlation of the entire distribution (doesn't make any sense upon review)
% Lesioned hemisphere
Coh_Ket1_Les = Coh_z_center_norm_sig(IX_Ket1_Les_M1,:);
Coh_Ket1_Les = Coh_Ket1_Les(:);
Coh_Ket2_Les = Coh_z_center_norm_sig(IX_Ket2_Les_M1,:);
Coh_Ket2_Les = Coh_Ket2_Les(:);
Coh_Ket3_Les = Coh_z_center_norm_sig(IX_Ket3_Les_M1,:);
Coh_Ket3_Les = Coh_Ket3_Les(:);
Coh_LDO1_Les = Coh_z_center_norm_sig(IX_LDO1_Les_M1,:);
Coh_LDO1_Les = Coh_LDO1_Les(:);

[rho(1),pval(1)] = corr(Coh_Ket1_Les,Coh_Ket2_Les,'Rows','pairwise');
[rho(2),pval(2)] = corr(Coh_Ket1_Les,Coh_Ket3_Les,'Rows','pairwise');
[rho(3),pval(3)] = corr(Coh_Ket2_Les,Coh_Ket3_Les,'Rows','pairwise');
[rho(4),pval(4)] = corr(Coh_LDO1_Les,Coh_Ket2_Les,'Rows','pairwise');

% plot scatter plots
figure; scatter(Coh_Ket1_Les,Coh_Ket2_Les)
lsline
hold on
figure; scatter(Coh_Ket1_Les,Coh_Ket3_Les)
lsline()
figure; scatter(Coh_Ket2_Les,Coh_LDO1_Les)
lsline()
title('Lesioned LDO&Ket; PreKetvsPostKet1(blue), PreKetvsPostKet2(red), Coh_z_sig_norm') 
% Unlesioned hemisphere
% Correlate pre and post ketamine periods of z scores
Coh_Ket1_UnLes_M1 = Coh_z_center_norm_sig(IX_Ket1_UnLes_M1,:);
Coh_Ket1_UnLes_M1 = Coh_Ket1_UnLes_M1(:);
Coh_Ket2_UnLes_M1 = Coh_z_center_norm_sig(IX_Ket2_UnLes_M1,:);
Coh_Ket2_UnLes_M1 = Coh_Ket2_UnLes_M1(:);
Coh_Ket3_UnLes_M1 = Coh_z_center_norm_sig(IX_Ket3_UnLes_M1,:);
Coh_Ket3_UnLes_M1 = Coh_Ket3_UnLes_M1(:);


[rho(1),pval(1)] = corr(Coh_Ket1_UnLes_M1,Coh_Ket2_UnLes_M1,'Rows','pairwise');
[rho(2),pval(2)] = corr(Coh_Ket1_UnLes_M1,Coh_Ket3_UnLes_M1,'Rows','pairwise');
[rho(3),pval(3)] = corr(Coh_Ket2_UnLes_M1,Coh_Ket3_UnLes_M1,'Rows','pairwise');

% plot scatter plots
figure; scatter(Coh_Ket1_UnLes_M1,Coh_Ket2_UnLes_M1)
lsline
hold on
scatter(Coh_Ket1_UnLes_M1,Coh_Ket3_UnLes_M1)
lsline()
title('UnLesioned LDO&Ket; PreKetvsPostKet1(blue), PreKetvsPostKet2(red), Coh_z_sig_norm') 

%%% Lesioned hemisphere %%%%%
% LDOPA
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_LDO1_Les_M1,:))) % Ldopa 1
title('Coherance z normalized Pre LDOPA Lesion M1')
colorbar
[out, vL1L, six] = sort_matrix(Coh_z_center_norm_sig(IX_LDO1_Les_M1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_LDO2_Les_M1,:))) % Ldopa 2
title('Coherance z normalized Post LDOPA Lesion M1')
colorbar
[out, vL2L, six] = sort_matrix(Coh_z_center_norm_sig(IX_LDO2_Les_M1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_LDO3_Les_M1,:))) % Ldopa 3
title('Coherance z normalized Peak 80 Lesion M1')
colorbar
[out, vL3L, six] = sort_matrix(Coh_z_center_norm_sig(IX_LDO3_Les_M1,:));

pubify_figure_axis
ax=gca;
ax.FontSize=35;

% Ketamine
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket1_Les_M1,:))) % Ket 1
title('Coherance z normalized Pre Ket Lesion M1')
colorbar
[out, vK1L, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket1_Les_M1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket2_Les_M1,:),'kmeans',7)) % Ket 2
title('Coherance z normalized kmeans Pre Ket Lesion M1')
colorbar
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket2_Les_M1,:))) % Ket 2
title('Coherance z normalized Post Ket Lesion M1')
colorbar
[out, vK2L, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket2_Les_M1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket3_Les_M1,:))) % Ket 3
title('Coherance z normalized Post Ket Late Lesion M1')
colorbar
[out, vK3L, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket2_Les_M1,:));

%%%%% Unlesioned hemisphere %%%%%
% LDOPA 
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_LDO1_UnLes_M1,:))) % LDOPA 1
title('Coherance z normalized Pre LDOPA UnLesion M1')
colorbar
[out, vL1UL, six] = sort_matrix(Coh_z_center_norm_sig(IX_LDO1_UnLes_M1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_LDO2_UnLes_M1,:),'kmeans',7)) % LDOPA 2
title('Coherance z normalized kmeans Pre LDOPA UnLesion M1')
colorbar
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_LDO2_UnLes_M1,:))) % LDOPA 2
title('Coherance z normalized Post LDOPA UnLesion M1')
colorbar
[out, vL2UL, six] = sort_matrix(Coh_z_center_norm_sig(IX_LDO2_UnLes_M1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_LDO3_UnLes_M1,:))) % LDOPA 3
title('Coherance z normalized Peak 80 UnLesion M1')
colorbar
[out, vL3UL, six] = sort_matrix(Coh_z_center_norm_sig(IX_LDO3_UnLes_M1,:));

% Ketamine
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket1_UnLes_M1,:))) % Ket 1
title('Coherance z normalized Pre Ket UnLesion M1')
colorbar
[out, vK1UL, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket1_UnLes_M1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket2_UnLes_M1,:),'kmeans',7)) % Ket 2
title('Coherance z normalized kmeans Pre Ket UnLesion M1')
colorbar
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket2_UnLes_M1,:))) % Ket 2
title('Coherance z normalized Post Ket UnLesion M1')
colorbar
[out, vK2UL, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket2_UnLes_M1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket3_UnLes_M1,:))) % Ket 3
title('Coherance z normalized Post Ket Late UnLesion M1')
colorbar
[out, vK3UL, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket3_UnLes_M1,:));

% plotting polar histogram for neuron 40 in Ket 2 group Lesion
rad_binsize = deg2rad(15);
rad_edges = -pi:rad_binsize:pi;
K2_Les = TBL(IX_Ket2_Les_M1,:);
d_Ket2_Les_40 = table2cell(K2_Les(40,18));
d_Ket2_Les_40 = cell2mat(d_Ket2_Les_40);
Ket2_Les_40 = table2cell(K2_Les(40,14));
Ket2_Les_40 = cell2mat(Ket2_Les_40);
clrs = lines(3);
figure;polarhistogram('BinEdges',rad_edges,'BinCounts',d_Ket2_Les_40(3,:) ...
    ,'DisplayStyle','stairs','EdgeColor','k')
hold on
polarhistogram('BinEdges',rad_edges,'BinCounts',Ket2_Les_40(3,:) ,'FaceAlpha',.2,'EdgeColor',clrs(1,:),'FaceColor',clrs(1,:))
Ket2_Les_40_peak = table2cell(K2_Les(40,34));                       
Ket2_Les_40_peak = cell2mat(Ket2_Les_40_peak);
peak_fq = find(Ket2_Les_40_peak(:,1)==max(Ket2_Les_40_peak)); % peak is at 74.5 Hz
% Plotting for gamma_80
figure;polarhistogram('BinEdges',rad_edges,'BinCounts',d_Ket2_Les_40(4,:) ...
    ,'DisplayStyle','stairs','EdgeColor','k')
hold on
polarhistogram('BinEdges',rad_edges,'BinCounts',Ket2_Les_40(4,:) ,'FaceAlpha',.2,'EdgeColor',clrs(1,:),'FaceColor',clrs(1,:))
% filtering LFp for 65 - 85 Hz range
% load('amp-A-056_LFP.mat')
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();
LFP = LK_Load_and_Clean_LFP('C:\Users\Stephen Cowen\Box\Cowen Laboratory\Data\LID_Ketamine_Single_Unit_R56\Rat320\20\LFP','amp-A-056_LFP.mat');
F = SPEC_create_filters([65 85],LFP.sFreq);
% Run the interval times in Q12
Ket_2.start = e_start_end_uS(1,1) + intervals_around_drug_min{1}(2,1)*60*1e6;
Ket_2.stop = e_start_end_uS(1,1) + intervals_around_drug_min{1}(2,2)*60*1e6;
ix1 = binsearch(LFP.t_uS , Ket_2.start);
ix2 = binsearch(LFP.t_uS , Ket_2.stop);
L = filtfilt(F{1},LFP.LFP(ix1:ix2));
TSr = Restrict(TS(25,1),LFP.t_uS(ix1),LFP.t_uS(ix2));
[PhC,phinfo] = SPEC_spike_field_coupling(TS,[LFP.t_uS(ix1:ix2), L],75,'thresh_prctile_of_power',80);
figure;polarhistogram('BinEdges',rad_edges,'BinCounts',PhC(25).sh_hist_rad_mn ...
    ,'DisplayStyle','stairs','EdgeColor','k')
hold on
polarhistogram('BinEdges',rad_edges,'BinCounts',PhC(25).hist_rad ,'FaceAlpha',.2,'EdgeColor',clrs(1,:),'FaceColor',clrs(1,:))
%%
%%% Control %%%%
% Ketamine
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket1_CM1,:))) % Ket 1
title('Coherance z normalized Pre Ket CM1')
colorbar
[out, vK1C, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket1_CM1,:));
% figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket2_CM1,:),'kmeans',7)) % Ket 2
% title('Coherance z normalized kmeans Pre CM1')
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket2_CM1,:))) % Ket 2
title('Coherance z normalized Post Ket CM1')
colorbar
[out, vK2C, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket2_CM1,:));
figure;imagesc(TBL.Coh_psd_fqs{1,1},[],sort_matrix(Coh_z_center_norm_sig(IX_Ket3_CM1,:))) % Ket 3
title('Coherance z normalized Post Ket Late CM1')
colorbar
[out, vK3C, six] = sort_matrix(Coh_z_center_norm_sig(IX_Ket3_CM1,:));

pubify_figure_axis
ax=gca;
ax.FontSize=35;

% Correlate pre and post ketamine periods of z scores
Coh_Ket1_CM1 = Coh_z_center_norm_sig(IX_Ket1_CM1,:);
Coh_Ket1_CM1 = Coh_Ket1_CM1(:);
Coh_Ket2_CM1 = Coh_z_center_norm_sig(IX_Ket2_CM1,:);
Coh_Ket2_CM1 = Coh_Ket2_CM1(:);
Coh_Ket3_CM1 = Coh_z_center_norm_sig(IX_Ket3_CM1,:);
Coh_Ket3_CM1 = Coh_Ket3_CM1(:);


[rho(1),pval(1)] = corr(Coh_Ket1_CM1,Coh_Ket2_CM1,'Rows','pairwise');
[rho(2),pval(2)] = corr(Coh_Ket1_CM1,Coh_Ket3_CM1,'Rows','pairwise');
[rho(3),pval(3)] = corr(Coh_Ket2_CM1,Coh_Ket3_CM1,'Rows','pairwise');

% plot scatter plots
figure; scatter(Coh_Ket1_CM1,Coh_Ket2_CM1)
lsline
hold on
scatter(Coh_Ket1_CM1,Coh_Ket3_CM1)
lsline()
title('Control Sal&Ket; PreKetvsPostKet1(blue), PreKetvsPostKet2(red), Coh_z_sig_norm') 
% fill in the psd fqs
fill_fqs = TBL.Coh_psd_fqs{1,1};
for ii = 1:Rows(TBL)
    if isempty(TBL.Coh_psd_fqs{ii,1})
       TBL.Coh_psd_fqs{ii,1} = single(fill_fqs);
    end
end

%% Plotting the ang p and ang z of all intervals
%%%% Ketamine interval 1 %%%%
% Lesion
M = TBL.Ang_p(IX_Ket1_Les_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket1_Les_M1,:);
MZ = TBL.Ang_z(IX_Ket1_Les_M1,:);
figure
subplot(141)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Peak 80/Pre Ket Lesion M1')
subplot(142)
Error_bars(MZ)
title('phase ang_z; Peak 80/Pre Ket Lesion M1')
% ttest
[h,p(1)] = ttest(MZ(:,3)-MZ(:,2)); % Comparing gamma to beta
[h,p(2)] = ttest(MZ(:,3)-MZ(:,4)); % comparing gamma to 80 Hz

bonf_holm(p)
% Unlesion
M = TBL.Ang_p(IX_Ket1_UnLes_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket1_UnLes_M1,:);
MZ = TBL.Ang_z(IX_Ket1_UnLes_M1,:);
subplot(143)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Peak 80/Pre Ket UnLesion M1')
subplot(144)
Error_bars(MZ)
title('phase ang_z; Peak 80/Pre Ket UnLesion M1')

%%%% Ketamine interval 2 %%%%
% Lesion
M = TBL.Ang_p(IX_Ket2_Les_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket2_Les_M1,:);
MZ = TBL.Ang_z(IX_Ket2_Les_M1,:);
figure
subplot(141)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Post Ket Lesion M1')
subplot(142)
Error_bars(MZ)
title('phase ang_z; Post Ket Lesion M1')
% Unlesion
M = TBL.Ang_p(IX_Ket2_UnLes_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket2_UnLes_M1,:);
MZ = TBL.Ang_z(IX_Ket2_UnLes_M1,:);
subplot(143)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Post Ket UnLesion M1')
subplot(144)
Error_bars(MZ)
title('phase ang_z; Post Ket UnLesion M1')

%%%% Ketamine interval 3 %%%%
% Lesion
M = TBL.Ang_p(IX_Ket3_Les_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket3_Les_M1,:);
MZ = TBL.Ang_z(IX_Ket3_Les_M1,:);
figure
subplot(141)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Post Ket Late Lesion M1')
subplot(142)
Error_bars(MZ)
title('phase ang_z; Post Ket Late Lesion M1')
% Unlesion
M = TBL.Ang_p(IX_Ket3_UnLes_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket3_UnLes_M1,:);
MZ = TBL.Ang_z(IX_Ket3_UnLes_M1,:);
subplot(143)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Post Ket Late UnLesion M1')
subplot(144)
Error_bars(MZ)
title('phase ang_z; Post Ket Late UnLesion M1')
%%
%%%% LDOPA interval 1 %%%%
% Lesion
M = TBL.Ang_p(IX_LDO1_Les_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_LDO1_Les_M1,:);
MZ = TBL.Ang_z(IX_LDO1_Les_M1,:);
figure
subplot(141)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Pre LDOPA Lesion M1')
subplot(142)
Error_bars(MZ)
title('phase ang_z; Pre LDOPA Lesion M1')
% Unlesion
M = TBL.Ang_p(IX_LDO1_UnLes_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_LDO1_UnLes_M1,:);
MZ = TBL.Ang_z(IX_LDO1_UnLes_M1,:);
subplot(143)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Pre LDOPA UnLesion M1')
subplot(144)
Error_bars(MZ)
title('phase ang_z; Pre LDOPA UnLesion M1')

%%%% LDOPA interval 2 %%%%
% Lesion
M = TBL.Ang_p(IX_LDO2_Les_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_LDO2_Les_M1,:);
MZ = TBL.Ang_z(IX_LDO2_Les_M1,:);
figure
subplot(141)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Post LDOPA Lesion M1')
subplot(142)
Error_bars(MZ)
title('phase ang_z; Post LDOPA Lesion M1')
% Unlesion
M = TBL.Ang_p(IX_LDO2_UnLes_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_LDO2_UnLes_M1,:);
MZ = TBL.Ang_z(IX_LDO2_UnLes_M1,:);
subplot(143)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Post LDOPA UnLesion M1')
subplot(144)
Error_bars(MZ)
title('phase ang_z; Post LDOPA UnLesion M1')

%%%% LDOPA interval 3 %%%%
% Lesion
M = TBL.Ang_p(IX_LDO3_Les_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_LDO3_Les_M1,:);
MZ = TBL.Ang_z(IX_LDO3_Les_M1,:);
figure
subplot(141)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Peak 80 Hz Lesion M1')
subplot(142)
Error_bars(MZ)
title('phase ang_z; Preak 80 Hz Lesion M1')
% Unlesion
M = TBL.Ang_p(IX_LDO3_UnLes_M1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_LDO3_UnLes_M1,:);
MZ = TBL.Ang_z(IX_LDO3_UnLes_M1,:);
subplot(143)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Peak 80 Hz UnLesion M1')
subplot(144)
Error_bars(MZ)
title('phase ang_z; Peak 80 Hz UnLesion M1')

%% Control Ketamine intervals
%%%% Ketamine interval 3 %%%%
M = TBL.Ang_p(IX_Ket1_CM1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket1_CM1,:);
MZ = TBL.Ang_z(IX_Ket1_CM1,:);
figure
subplot(121)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Pre Ket Control M1')
subplot(122)
Error_bars(MZ)
title('phase ang_z; Pre Ket Control M1')

%%%% Ketamine interval 2 %%%%
M = TBL.Ang_p(IX_Ket2_CM1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket2_CM1,:);
MZ = TBL.Ang_z(IX_Ket2_CM1,:);
figure
subplot(121)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Post Ket Control M1')
subplot(122)
Error_bars(MZ)
title('phase ang_z; Post Ket Control M1')

%%%% Ketamine interval 3 %%%%
M = TBL.Ang_p(IX_Ket3_CM1,:);
M_shuf = TBL.Ang_to_shuf_p(IX_Ket3_CM1,:);
MZ = TBL.Ang_z(IX_Ket3_CM1,:);
figure
subplot(121)
bar(nansum(M<.05)-nansum(M_shuf<.05))
title('<.05 p value neurons (shuffle subtracted) of phase ang; Post Ket Late Control M1')
subplot(122)
Error_bars(MZ)
title('phase ang_z; Post Ket Late Control M1')

%% Z score normalised Ang_z
M = TBL.Ang_z(IX_Ket_2,:);
M_n = Z_scores(M')';
figure;imagesc(M_n)
figure;Error_bars(M_n)

% Is the 50 Hz during ketamine different from teh other bands
% Within subjects paired ttest for the raliegh z scores for the diff
% between 50-8 and 50-80. Correct for bonferroni
[h,p(1)] = ttest(M(:,3)-M(:,2));
[h,p(2)] = ttest(M(:,3)-M(:,4));

bonf_holm(p)



figure;imagesc(TBL.Coh_circ_rtest_z)

%% Plotting for the spike field coupling for peak fqs code
p_th = 0.05;
for iN = 1:length(PhPk)
        if PhPk(iN).circ_p_to_shuff < p_th
            figure
            subplot(1,2,1)
            scatter(PhPk(iN).pk_freq,PhPk(iN).pk_pow)
            lsline
            ylabel('pow')
            xlabel('freq');
            axis tight
            title(sprintf('N%d, pts %1.3f ,  p %1.3f, z %1.3f',iN, PhPk(iN).circ_p_to_shuff ,PhPk(iN).circ_rtest_p, PhPk(iN).circ_rtest_z));
            subplot(1,2,2)
            polarhistogram(PhPk(iN).pk_phase,40)
            
        end
end