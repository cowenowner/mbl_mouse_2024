fqs = 1:.5:170;
E = LK_Load_Events_From_Excel('Event_times.xlsx');

if any(E.EventID=='LDOPAInjectionStart')
    t_start_end_min(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')+38;
    t_start_end_min(2) = t_start_end_min(1) + 20;
    sig_name = 'gamma_80';
elseif  any(E.EventID=='KetInjectionStart')
    t_start_end_min(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')+2;
    t_start_end_min(2) = t_start_end_min(1) + 20;
    sig_name = 'gamma_50';
end

t_start_end_uS = t_start_end_min*60e6;


lfp_dir = 'H:\LFP_analysis\Rat380\11';
[~,LFP_files] = find_files(fullfile(lfp_dir,'amp*.mat'));

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

RR_left = L(:,2) - L(:,1);
RR_right = L(:,5) - L(:,4);

figure
subplot(131)
pwelch(L(:,2),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
sgtitle('Best channel signal left hem ketamine')
subplot(132)
pwelch(L(:,1),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
sgtitle('700uM away signal left hem ketamine')
subplot(133)
pwelch(RR_left,LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
sgtitle('re reffed signal left hem ketamine')

figure
subplot(131)
pwelch(L(:,5),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
sgtitle('Best channel signal right hem ketamine')
subplot(132)
pwelch(L(:,4),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
sgtitle('700uM away signal right hem ketamine')
subplot(133)
pwelch(RR_right,LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
sgtitle('re reffed signal right hem ketamine')

figure
subplot(121)
pwelch(L(:,3),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
title('Best channel signal left hem Striatum ketamine')
subplot(122)
pwelch(L(:,6),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
title('Best channel signal right hem Striatum ketamine')

figure
subplot(121)
pwelch(L(:,3),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
title('Best channel signal left hem Striatum peak 80')
subplot(122)
pwelch(L(:,6),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
title('Best channel signal right hem Striatum peak 80')

%%
spec_win_sec = 5;
[S,w,t] = spectrogram(LFP.LFP,LFP.sFreq*spec_win_sec,LFP.sFreq*(spec_win_sec/2),fqs, LFP.sFreq);
    sFreq_of_spectrogram = 1/median(diff(t));
    t_min = t/60;
    S = 10*log10(abs(S))';
    %     p = prctile(S,[.5 99.5]);
    %     BIX = S<p(1,:);
    %     BIX = BIX | S>p(2,:);
    %     S(BIX) = nan;
%     Sn = movmedian(S,30,'omitnan');
    han_size_sec = 40;
    han_size_points = han_size_sec*sFreq_of_spectrogram;
    Sn = conv_filter(S,hanning(han_size_points)/sum(hanning(han_size_points)));
    figure
    imagesc(t_min,fqs,Sn')
    axis xy
    ylabel('Hz')
    colorbar
    caxis(prctile(Sn(:),[1 98.6]))
    
%%
figure
for ii = 1:Cols(L)
    subplot(ceil(Cols(L)/8),8,ii)
    pwelch(L(:,ii),LFP.sFreq,LFP.sFreq/2,256,LFP.sFreq)
    title(LFP_files{ii},'FontSize',7)
end