function [cc,lags,ccz] = xcorr_cowen(v1, v2, n_lags, varargin)
% this variant of xcorr is meant to overcome some of the problems and
% assumptions of the built in xcorr. 
% subtracts mean from v1 and v2
% provides a shuffle mean and std estimate for evaluating your cc
%
% See Cross_corr_demo
% For a list of some of the corruptions of cross-corrs and the solutions
% (via the generation of null distributions via suffling and randomizing in
% clever ways) see... (these are for spike trains, but the general
% principles and confounds are mostly the same for cross corss between
% continuous signals).
%
% Brody, C.D., 1999. Correlations without synchrony. Neural Comput 11, 1537–1551. https://doi.org/10.1162/089976699300016133
% Palm, G., Aertsen, A.M., Gerstein, G.L., 1988. On the significance of correlations among neuronal spike trains. Biol Cybern 59, 1–11. https://doi.org/10.1007/BF00336885
%
%%%%%%%
% Cowen
n_shuffle = 500;
type = 'coeff'; % 'none' (default) | 'biased' | 'unbiased' | 'normalized' | 'coeff'
shuff_method = 'circshift'; % perm
 % shuff_method = 'permute'; % perm

Extract_varargin;

%% %%%%%%%%%%%%%%%%%%%%%%%
GIX = ~isnan(v1) |~isnan(v2) ;
v1 = v1(GIX) - mean(v1(GIX));
v2 = v2(GIX) - mean(v2(GIX));
v1 = v1(:); v2 = v2(:);
[cc,lags] = xcorr(v1 - mean(v1), v2 - mean(v2), n_lags, type);
cc = cc';
shuffle_xc = nan(n_shuffle,2*n_lags+1);
for iShuff = 1:n_shuffle
    switch shuff_method
        case 'permute'
            v1sh = v1(randperm(length(v1)));
        case 'circshift'
            v1sh = circshift(v1,round((rand(1,1)-.5)*.25*length(v1)));
        otherwise
            error('wtf')
    end
    % We could shuffle v2 as well but that is probably not necessary.
    shuffle_xc(iShuff,:) = xcorr(v1sh - mean(v1sh), v2 - mean(v2), n_lags, type)';
end

ccz = (cc-mean(shuffle_xc))./std(shuffle_xc);
%%
if nargout == 0
    figure
    subplot(1,2,1)
    plot(lags,cc);hold on
    plot(lags,mean(shuffle_xc),'r')
    plot(lags,mean(shuffle_xc) + 1.96*std(shuffle_xc),'r:')
    plot(lags,mean(shuffle_xc) - 1.96*std(shuffle_xc),'r:')
    ylabel(type)
    axis tight
    plot_vert_line_at_zero
    title('cc with 95% shuf CI')

    subplot(1,2,2)
    plot(lags,ccz)
    ylabel('z')
    axis tight
    plot_vert_line_at_zero
end