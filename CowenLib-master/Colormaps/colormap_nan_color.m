function colormap_nan_color(Data_Array,nan_color, cmap)
% adjusts the colormap so that nans now have their own color.
% data array is the original data used in teh imagesc plot.
% HOW TO CALL: first do imagesc, then call this function to adjust teh
% colormap.
%
% Cowen 2020. Found on matlab discussion forum - very clever.
if nargin < 2
    nan_color = [1 1 1];
end

if nargin >= 3
    colormap(cmap)
end

imAlpha=ones(size(Data_Array));
imAlpha(isnan(Data_Array))=0;
imagesc(Data_Array,'AlphaData',imAlpha);

set(gca,'color',nan_color);
