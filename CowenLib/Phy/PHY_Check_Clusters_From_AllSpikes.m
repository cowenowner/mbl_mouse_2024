function PHY_Check_Clusters_From_AllSpikes(SP, save_png)
% PHY_Check_Clusters_From_AllSpikes
% Cowen 2024
if nargin < 1
    load('AllSpikes.mat')
end
if nargin < 2
    save_png = false;
end

for iF = 1:length(SP)
    if any(diff(SP(iF).t_uS)<=0)
        error('There are some corrupted timestamps.')
    end
    figure(101010)
    clf
    subplot (3,3,1)
    % PLOT a subset of the 3 max amp waveforms.
    mx = max(SP(iF).WV.mWV);
    [~,six] = sort(mx,'descend');
    for ii = 1:3
        plot(SP(iF).WV.mWV(:,six(ii))+80*(ii-1),'-')
        hold on
    end
    title(sprintf('ID: %d, %s', SP(iF).cluster_id, SP(iF).PHYLabel))
    axis tight
    xlabel('sample')
    box off

    subplot (3,3,2)
    template_v = sum(abs(SP(iF).template),2);
    [~,six] = sort(template_v,'descend');
    plot(SP(iF).template(six(1:6),:)'); 
    axis tight
    title('1st 6 templts')
    xlabel('sample')
    axis off
  
    subplot (3,3,3)
    [b,x] = AutoCorr(SP(iF).t_uS/100,4,100);
    plot(x,b)
    xlabel('ms')
    axis tight
    box off

    title(sprintf('AC 4ms bin, clu %d nSpk %d',SP(iF).cluster_id,SP(iF).n_spikes))

    subplot (3,3,4)
    HistISI(SP(iF).t_uS/100)
    box off

    title(sprintf('Dur = %2.2f hr, %d/%d < 2ms',(SP(iF).t_uS(end) - SP(iF).t_uS(1))/3600e6,sum(diff(SP(iF).t_uS/1000) < 2), length(SP(iF).t_uS)))
    subplot(3,3,5)
    if 0
        PKS = SP(iF).template_features;
        dims = [2 1;2 3;4 3;4 1];
        signs = [1 1; 1 -1; -1 -1; -1 1];
        for iD = 1:Rows(dims)
            d = dims(iD,:);
            s = signs(iD,:);
            plot(PKS(1:5:end,d(1))*s(1),PKS(1:5:end,d(2))*s(2),'k.','MarkerSize',1)
            hold on
            plot(PKS_cell(:,d(1))*s(1),PKS_cell(:,d(2))*s(2),'r.')
        end
        xmin = pks_rng(2,4)*-1;
        xmax = pks_rng(2,2);
        ymin = pks_rng(2,3)*-1;
        ymax = pks_rng(2,1);
        set(gca,'XLim',[xmin xmax])
        set(gca,'YLim',[ymin ymax])
        plot_vert_line_at_zero
        plot_horiz_line_at_zero
    end

    if isfield(SP,'template_features')

        plot(SP(iF).template_features(:,1),SP(iF).template_features(:,2),'.')

    end
    title(sprintf('Depth %1.1f uM',SP(iF).neuropixels_depth_uM))

    subplot(3,3,6)
    [b,x] = AutoCorr(SP(iF).t_uS/100,20,100);
    plot(x,b)
    xlabel('ms')
    axis tight
    box off

    title('AC 20ms bin, 2sec')


    subplot (3,3,7:9)
    plot(SP(iF).t_uS/3600e6,SP(iF).amp_all,'.')
    hold on
    xlabel('hours')
    ylabel('amp')
    axis tight
    set(gcf,'Position',[554.6000  244.2000  812.8000  680.0000])
    if save_png
         saveas(gcf,fullfile(pwd,[nm '.png']))
    end
end
