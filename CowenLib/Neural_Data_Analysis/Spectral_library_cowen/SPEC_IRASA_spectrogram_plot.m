function [Spec] = SPEC_IRASA_spectrogram_plot(Spec)
% FROM: https://purr.purdue.edu/publications/1987/1
% Also see https://pubmed.ncbi.nlm.nih.gov/36639900/
figure
subplot(4,1,1);
imagesc([],Spec.freq,log(Spec.frac))
colorbar
title(sprintf('log Spec.frac B=%1.2f C=%1.2f ', mean(Spec.Beta),mean(Spec.Cons)))

subplot(4,1,2)
plot(Spec.Beta); hold on
plot(Spec.Cons)
colorbar
legend('Beta','Cons'); legend boxoff
axis tight

subplot(4,1,3);
imagesc([],Spec.freq,Spec.osci) % is it abs() or is it real()?????? In the Example.m, it does NEITHER!
colorbar
title('abs Spec.osci')

% subplot(3,1,3);
% loglog(Spec.freq,Spec.mixd,'b'); hold on;
% loglog(Spec.freq,Spec.frac,'r');
% legend('mixd','frac');
% 
% subplot(3,1,4);
% loglog(Spec.freq,mean(Spec.mixd,2),'b'); hold on
% loglog(Spec.freq,mean(Spec.frac,2),'r');
% legend('mixd','frac');

subplot(4,1,4);
plot(Spec.freq, mean(Spec.osci,2));
hold on
plot(Spec.freq, mean(Spec.frac,2));
axis tight
title('mean(Frac.osci,2) mean(Spec.frac,2)')
plot_horiz_line_at_zero;
legend('osci','frac');

%% pc exploratory.
[pc,sc] = pca([Spec.osci' Spec.frac' Spec.Beta Spec.Cons]);
figure
subplot(2,2,1)
plot(sc(:,1),sc(:,2),'.')
subplot(2,2,2)
plot(sc(:,3),sc(:,4),'.')
subplot(2,2,3)
plot(sc(:,5),sc(:,6),'.')
subplot(2,2,4)
plot(sc(:,7),sc(:,8),'.')
