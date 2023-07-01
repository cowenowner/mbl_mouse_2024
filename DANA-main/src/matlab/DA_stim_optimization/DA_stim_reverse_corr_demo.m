% A demonstration of how we could, use reverse correlations to determine
% how sequences of stims affect dopamine release.
% Two approaches here. In one, we just look at say to sequential ISIs. Vary these 2
% ISIs randomly throuhg the course of the experiment and then map the
% sequence to dopamine output.
% In the second, we give randome pulses say in groups of 5 and do this
% again and again and interleave yoked fixed frequency stimulation. 
% Then split the trials into those with say low, medium,
% and high dopamine release. Once grouped, look at the mean ISI's that lead
% up to high vs. low release and also look to see how regular those ISIs
% were and how they compared to fixed-frequencies of the same duration.

%
%% Cowen 2020
%% Each element of TGT is the delay, in sec, between adjacten stims (sec).
close all
mean_freq = 20; % mean stim frequency.
n_ISIs_in_tgt = 2; % 3 pulses for now as hard to visualise > 3 for reverse correlations.
n_ISIs_total = 100; % 3 pulses for now as hard to visualise > 3 for reverse correlations.
rng(4); % For replication.
rand_factor = .001;
TGT = []; % the targets that will produce the most dopamine
% By having >1 TGT, you make the problem much harder, but probably more
% realistic. If it does fluctuate, then how do we get it to better converge
% on one solution instead of some half-assed middle solution?
%
TGT(1,:) = rand(1,n_ISIs_in_tgt) /mean_freq; % Random sequences
TGT(1,:) = [0.018351     0.078274]; % Random sequences
STIM = rand(1,n_ISIs_total+n_ISIs_in_tgt)*1/mean_freq;
for ii = 1:length(STIM)-n_ISIs_in_tgt
    v = TGT*STIM(ii:ii + n_ISIs_in_tgt - 1)';
    DA(ii) = -1*sum((v).^1.2) + rand(1,1)*rand_factor;
%     DA(ii) = -1*sum(sin(v*50).^2) + rand(1,1)*rand_factor;
end
% Some sort of meshgrid thing to interp dopamine to different stim levels.
x = STIM(1:end-1);
y = STIM(2:end);
nrec = length(x)-2;
x = x(1:nrec)'; y = y(1:nrec)'; DA = DA(1:nrec)';
% Generate a fit and surface plot using a higher order model.
%sf = fit([x(1:nrec)',y(1:nrec)'],DA(1:nrec)','poly23');
xlin=linspace(min(x),max(x),33)';        % Create x,y linear space
ylin=linspace(min(y),max(y),33)';
[X,Y]=meshgrid(xlin,ylin);              % Create mesh [x y]
Z = griddata(x,y,DA,X,Y,'cubic');          % Interpolate with bicubic functions            
%
figure
subplot(1,2,1)
mesh(X,Y,Z); % interpolated             % Fancy plots for demosntration
hold on
plot3(x,y,DA,'.','MarkerSize',15)
xlabel('ISI1 sec'),ylabel('ISI2 sec'),zlabel('DA')

subplot(1,2,2)
imagesc(xlin, ylin, Z)
xlabel('ISI1 sec'),ylabel('ISI2 sec')
axis xy
colorbar
axis square

% STA

