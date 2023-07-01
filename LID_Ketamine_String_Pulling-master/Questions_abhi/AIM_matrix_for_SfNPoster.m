
%%
% Matrices of AIM scores for all LDOPA sessions for Rat 320
% Rat320_Ses3_L = [ 0 0; 20 .38; 40 1.5; 59 1.75; 80 2.5; 100 2.63; 120 1.75; 140 1; 160 .25; 180 0];
% Rat320_Ses8_L = [ 0 0; 20 1; 40 2.5; 59 2.63; 80 2; 100 3.25; 120 2.25; 140 1.25; 160 .38; 180 0];
% Rat320_Ses9_L = [ 0 0; 20 1; 40 1.63; 59 3; 80 2.38; 100 3.25; 120 2; 140 1; 160 .25; 180 0];
% Rat320_Ses13_L = [ 0 0; 20 1.88; 40 2.75; 59 3.5; 80 3; 100 2.25; 120 1.75; 140 1.63; 160 .25; 180 .25];
% Rat320_Ses15_L = [ 0 0; 20 1.25; 40 2.25; 59 3; 80 2.13; 100 2.25; 120 2.25; 140 1.5; 160 .25; 180 0];
% Rat320_Ses18_L = [ 0 0; 20 1; 40 2.75; 59 1.75; 80 3; 100 2; 120 1.25; 140 .25; 160 0; 180 0.5];
% Rat320_Ses21_L = [ 0 0; 20 1.25; 40 2.75; 59 2; 80 3; 100 1.88; 120 1.75; 140 1.5; 160 .5; 180 0.25];
Rat320_Ses3_L = [ -60 0; -40 .33; -20 1.67; -1 1.67; 20 3; 40 3.33; 60 2; 80 1; 100 .33; 120 0];
Rat320_Ses8_L = [ -60 0; -40 1; -20 2; -1 2.17; 20 2; 40 3.67; 60 2.33; 80 1.33; 100 .5; 120 0];
Rat320_Ses9_L = [ -60 0; -40 1; -20 1.83; -1 3.33; 20 2.5; 40 3.67; 60 2; 80 1.33; 100 .33; 120 0];
Rat320_Ses13_L = [ -60 0; -40 1.5; -20 3; -1 3.33; 20 3.33; 40 2.33; 60 1.67; 80 2; 100 .33; 120 .33];
Rat320_Ses15_L = [ -60 0; -40 1; -20 1.67; -1 3.33; 20 2.17; 40 2; 60 2.33; 80 1.33; 100 .33; 120 0];
Rat320_Ses18_L = [ -60 0; -40 1.33; -20 2.67; -1 1; 20 3.33; 40 2.33; 60 1.67; 80 .33; 100 0; 120 0.67];
Rat320_Ses21_L = [ -60 0; -40 1; -20 3; -1 2.33; 20 3.67; 40 1.83; 60 2.33; 80 1.83; 100 .67; 120 0.33];
% Concatenating matrices and getting the mean and SEM for all timepoints
Cat_allses = cat(3,Rat320_Ses3_L,Rat320_Ses8_L,Rat320_Ses9_L,Rat320_Ses13_L,Rat320_Ses15_L, ...
   Rat320_Ses18_L,Rat320_Ses21_L);
Avg_allses = mean(Cat_allses,3);
SEM_allses = std(Cat_allses,[],3)./sqrt(size(Cat_allses,3));
% Std_allses = std(Cat_allses,[],3);
% plotting the avg and sem
figure
errorbar(Avg_allses(:,1),Avg_allses(:,2),SEM_allses(:,2))
hold on

Cat_allses = cat(2,Rat320_Ses3_L(:,2),Rat320_Ses8_L(:,2),Rat320_Ses9_L(:,2),Rat320_Ses13_L(:,2),Rat320_Ses15_L(:,2), ...
   Rat320_Ses18_L(:,2),Rat320_Ses21_L(:,2));
%% 
% Creating matrices for AIM scores for all LDOPA and Ket sessions for Rat
% 320
% Rat320_Ses4_LK = [ 0 0; 20 2.38; 40 1.5; 59 2.13; 80 1.66; 100 1.88; 120 1.38; 140 1.5; 160 .5; 180 0];
% Rat320_Ses6_LK = [ 0 0; 20 1.25; 40 1.88; 59 2; 80 1.75; 100 1.38; 120 1.75; 140 1.75; 160 .75; 180 0];
% Rat320_Ses11_LK = [ 0 0; 20 .63; 40 2.88; 59 2.75; 80 1; 100 1.63; 120 1.13; 140 1.88; 160 1.75; 180 .25];
% Rat320_Ses14_LK = [ 0 0; 20 1.38; 40 2.25; 59 3.13; 80 1; 100 .75; 120 1.38; 140 .63; 160 .5; 180 .25];
% Rat320_Ses17_LK = [ 0 0; 20 1.63; 40 2; 59 3.38; 80 1; 100 1.25; 120 1.13; 140 1; 160 .25; 180 0];
% Rat320_Ses20_LK = [ 0 0; 20 2.13; 40 2.5; 59 2.25; 80 1; 100 1.63; 120 1.5; 140 .75; 160 1; 180 0];
% Rat320_Ses23_LK = [ 0 0; 20 .5; 40 2.88; 59 3.25; 80 1.5; 100 1.25; 120 1.5; 140 1; 160 .5; 180 0];
% Rat320_Ses24_LK = [ 0 0; 20 .88; 40 1; 59 2.88; 80 .75; 100 1; 120 1.25; 140 1; 160 .63; 180 0];
Rat320_Ses4_LK = [ -60 0; -40 2.5; -20 1.67; -1 2.17; 20 1.5; 40 1.83; 60 1.67; 80 1.33; 100 .67; 120 0];
Rat320_Ses6_LK = [ -60 0; -40 1; -20 1.5; -1 1.33; 20 1; 40 1.5; 60 1.67; 80 2; 100 .67; 120 0];
Rat320_Ses11_LK = [ -60 0; -40 .83; -20 3.17; -1 3.33; 20 .67; 40 .83; 60 .83; 80 1.83; 100 1.67; 120 .33];
Rat320_Ses14_LK = [ -60 0; -40 1.67; -20 1.67; -1 3.17; 20 .67; 40 1; 60 1.67; 80 .83; 100 .67; 120 .33];
Rat320_Ses17_LK = [ -60 0; -40 1.83; -20 1.67; -1 3.17; 20 .67; 40 1; 60 1.33; 80 1.33; 100 .33; 120 0];
Rat320_Ses20_LK = [ -60 0; -40 2.17; -20 3; -1 2.67; 20 .67; 40 1.67; 60 1.33; 80 1; 100 1; 120 0];
Rat320_Ses23_LK = [ -60 0; -40 .67; -20 3.17; -1 3; 20 1.33; 40 1; 60 1.33; 80 1.33; 100 .67; 120 0];
Rat320_Ses24_LK = [ -60 0; -40 1; -20 .33; -1 3.17; 20 .67; 40 1; 60 1; 80 1.33; 100 .83; 120 0];
% Concatenating matrices and getting the mean and SEM for all timepoints
Cat_allses = cat(3,Rat320_Ses4_LK,Rat320_Ses6_LK,Rat320_Ses11_LK,Rat320_Ses14_LK,Rat320_Ses17_LK, ...
   Rat320_Ses20_LK,Rat320_Ses23_LK,Rat320_Ses24_LK);
Avg_allses_ket = mean(Cat_allses,3);
SEM_allses = std(Cat_allses,[],3)./sqrt(size(Cat_allses,3));
% Std_allses = std(Cat_allses,[],3);
% plotting the avg and sem
errorbar(Avg_allses_ket(:,1),Avg_allses_ket(:,2),SEM_allses(:,2))
pubify_figure_axis
plot_ref_line(0);

Cat_allses_ket = cat(2,Rat320_Ses4_LK(:,2),Rat320_Ses6_LK(:,2),Rat320_Ses11_LK(:,2),Rat320_Ses14_LK(:,2),Rat320_Ses17_LK(:,2), ...
   Rat320_Ses20_LK(:,2),Rat320_Ses23_LK(:,2),Rat320_Ses24_LK(:,2));

%% ttest for 20 and 40 min scores post ketamine
Ket_20 = [1.5 1 .67 .67 .67 .67 1.33 .67];
Sal_20 = [3 2 2.5 3.33 2.17 3.33 3.67];
trying = Sal_20 - Ket_20;
% ;1.22875000000000
% ;2.73714285714286]
[h,pval] = ttest2(Sal_20,Ket_20,'Vartype','unequal');
% pval = signrank(Avg_allses([5 6 7],2),Avg_allses_ket([5 6 7],2));
for iA = 1:length(Cat_allses)
    pval(iA) = ranksum(Cat_allses(iA,:),Cat_allses_ket(iA,:));
end


