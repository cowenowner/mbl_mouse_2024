function OUT = Corr_upper_only(Q)
% function OUT = Corr_upper_only(Q)
%
% Computes the correlation coefficient for pairs of cols in a matrix Q and
% only returns the upper diagonal since the data is redundant for the lower
% diagonal and the diagonal is 1.
%
% Cowen 2023
[R,P] = corr(Q);
ONES = ones(size(R));
UPPER = triu(ONES,1) > 0;
[i,j]=ind2sub(size(UPPER),find(UPPER));
OUT.r = R(UPPER);
OUT.p = P(UPPER);
OUT.i = i;
OUT.j = j;

UPPER = double(UPPER);
UPPER(UPPER==0) = nan;

OUT.Rfull = R.*UPPER;

if nargout == 0
    % Make a pretty plot
    figure
    subplot(1,2,1)
    histogram(OUT.r,100);
    pubify_figure_axis
    axis square
    xlabel('r')

    subplot(1,2,2)
    imagesc(OUT.Rfull);
    colorbar('eastoutside')
    colormap(turbo)
    axis square
    pubify_figure_axis
end