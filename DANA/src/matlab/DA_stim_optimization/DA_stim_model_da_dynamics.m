% Reframing the problem here from optimizing to modeling.
% 
% Assumptions: mPFC stim produces a very non-linear looking dopamine
% response in the NAc. In the first test, this response is a linear
% convolution the multi-modal DA release with the pulse sequance. Once we
% have a system that can model this, then we can construct more complex
% functions of input and output.
%
%%
x = 0:.01:.6;
dur_s = 60*1;
stim_mean_hz = 20;
pulse_dur_s = 0.005;

% DA_response = wblpdf(x,.5,1.5) + normpdf(x+1,.5,3)
% DA_response = wblpdf(x,.2,1.5) + wblpdf(x-.4,.2,1.1); % Create some srange arbitrary non-linear, non-monotonic dopamine response.
DA_response = wblpdf(x,.2,1.5) + normpdf(x,.4,.06); % Create some srange arbitrary non-linear, non-monotonic dopamine response.
% DA_response = wblpdf(x,.2,1.5) ; % Create some srange arbitrary non-linear, non-monotonic dopamine response.
DA_response = DA_response/sum(DA_response);
% pad with zeros in front so that the events sync.
% DA_response = [zeros(size(DA_response)) DA_response];
% plot(x, wblpdf(x,.5,1.5))
figure
plot(x,DA_response); 

ylabel('DA');xlabel('s')

DA_response_pad = [zeros(size(DA_response)) DA_response];

%% Create the spike train.
x_s = linspace(0,dur_s,50000);
npts_in_stim = find(x_s > pulse_dur_s,1,'first');
S = zeros(size(x_s));
ISIs = abs(normrnd(1/stim_mean_hz,1/stim_mean_hz,1,10000));
% ISIs(1) = 0;
t_s = cumsum(ISIs);
t_s = t_s(t_s<dur_s);
ix = dsearchn(x_s(:),t_s(:));
S(ix) = 1;
S = conv(S,ones(1,npts_in_stim),'same');
% DA = convn(S,DA_response(end:-1:1),'same');
DA = convn(S,DA_response_pad,'same');
% Add noise
DA = DA + randn(size(DA))*mean(DA)/2;
% Add non-stationarity.
x = linspace(0,2*pi*5,length(DA));
sinwave = (sin(x)+1)/7;
DA = DA + sinwave;

figure
plot(DA+sinwave)

figure
bar(x_s,S)
hold on
plot(x_s,DA)
%% Now that we have the train, can we figure out the kernel response?
TR_IX = x_s<dur_s/2;
TE_IX = x_s>dur_s/2;
%
numFeatures = 1;
numResponses = 1;
numHiddenUnits = length(DA_response)+10;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];

% maxEpochs = 60;
% miniBatchSize = 20;

options = trainingOptions('adam', ...
    'MaxEpochs',200, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');
% Detrend
% Shan = convn(S,hanning(length(DA_response)+210),'same');
% Sn = zscore(Shan);
% DAn = zscore(DA);

Sn = S;
DAn = DA;

net = trainNetwork(Sn(TR_IX),DAn(TR_IX),layers,options);
% net = predictAndUpdateState(net,XTrain);
% [net,YPred] = predictAndUpdateState(net,YTrain(end));

[DApred] = predict(net,Sn(TE_IX));

figure
plot(x_s(TR_IX),Sn(TR_IX),'b',x_s(TE_IX),Sn(TE_IX),'r')
hold on
plot(x_s(TR_IX),DAn(TR_IX),'g',x_s(TE_IX),DAn(TE_IX),'c')

figure
bar(x_s(TE_IX),Sn(TE_IX))
hold on
plot(x_s(TE_IX),DAn(TE_IX),'k')
plot(x_s(TE_IX),DApred,'r')


%%%%%%%%%%%%%%%%%%%%%%%%%
%% TDNN: use timedelaynet or narxnet
%%%%%%%%%%%%%%%%%%%%%%%%%
% inputDelays = 1:length(DA_response)+10;
%%
X = num2cell(Sn(TR_IX)); % seems dumb, but narx and timedelay seem to net need cell arrays.
Xtest = num2cell(Sn(TE_IX)); % seems dumb, but narx and timedelay seem to net need cell arrays.
T = num2cell(DAn(TR_IX));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputDelays = 1:2;
inputDelays = 1:length(DA_response)+10;
% hiddenSizes = 1:length(DA_response)+10;
hiddenSizes = 10;
tdnn_net = timedelaynet(inputDelays,hiddenSizes);
[Xs,Xi,Ai,Ts] = preparets(tdnn_net,X,T);
tdnn_net_trained = train(tdnn_net,Xs,Ts,Xi,Ai);
view(tdnn_net_trained)
out = tdnn_net_trained(Xtest);
out = cell2mat(out);

figure
bar(x_s(TE_IX),Sn(TE_IX))
hold on
plot(x_s(TE_IX),DAn(TE_IX),'k')
plot(x_s(TE_IX),out,'r')
