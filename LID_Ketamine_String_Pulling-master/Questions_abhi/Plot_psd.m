data_dir = 'E:\SfN_LFP';
sFreq = 500;
ana_range_min = [2 30];
ana_t_sec = ana_range_min(1)*60:1/sFreq:ana_range_min(2)*60;
conditions_to_compare = {'Sal_Ket_control_lesion_unlesion_m1' 'Sal_Ket_control_lesion_unlesion_striatum'};
for iD = 1:length(conditions_to_compare)
    files = dir(fullfile(data_dir,conditions_to_compare{iD},'*.mat'));
    for iF = 1:length(files)
        D = load(fullfile(data_dir,conditions_to_compare{iD},files(iF).name));
        t_sec = (0:(length(D.LFP.data)-1))/D.LFP.LFP_sFreqj;
        if any(contains(fieldnames(D.LFP),'Sal_Start_min'))
            t_sec = t_sec - D.LFP.Sal_Start_min*60;
        else
            t_sec = t_sec - D.LFP.Ket_Start_min*60;
        end
        LFP = single(interp1(t_sec,double(D.LFP.data)*D.LFP.to_uV_conversion,ana_t_sec));
        LFP = [ana_t_sec(:) LFP(:)];
        %plot
        figure;
        [pxx,f] = pwelch(LFP(:,2),sFreq,sFreq/2,1:.25:200,sFreq);
        plot(f,10*log10(pxx),'LineWidth',2.0)
        pubify_figure_axis
        set(gca,'ylim',[-10 40])
        title(files(iF).name,'FontSize',7)
        grid off
        
    end
end

pow_psd = 10*log10(pxx);