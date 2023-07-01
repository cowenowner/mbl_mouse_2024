

for iN = 1:size(PhC)
%       if PhC(iN).Ang_p(4) < 0.05 && PhC(iN).Ang_to_shuf_p(4) > 0.05
         figure
%          subplot(2,ceil((Cols(DATA)-1)/2),iF)
         polarhistogram('BinEdges',rad_edges,'BinCounts',PhC(iN).sh_hist_rad_mn(4,:) ,'DisplayStyle','stairs','EdgeColor','k')
         hold on
         polarhistogram('BinEdges',rad_edges,'BinCounts',PhC(iN).sh_hist_rad_mn(4,:)  + PhC(iN).sh_hist_rad_95ci(4,:)  ,'DisplayStyle','stairs','EdgeColor',[.5 .5 .5])
         polarhistogram('BinEdges',rad_edges,'BinCounts',PhC(iN).hist_rad(4,:) ,'FaceAlpha',.5)
         title(sprintf('Neuron %d p=%1.3f z=%2.1f pshuff=%1.3f ',iN,PhC(iN).Ang_p(4),PhC(iN).Ang_z(4),PhC(iN).Ang_to_shuf_p(4)),'FontSize',10)
%      end
 end