function firings = Izhikevich_model_neurons(secs_of_data, Ne, Ni, neuron_params, fq)
% Created by Eugene M. Izhikevich, February 25, 2003
% modified heavily by cowen 2016
%
% OUTPUT: first column is the time (ms), the second is cell ID.
%
% Izhikevich_model_neurons(60,40,20,'low_firing_rate',[2,4;5,2]);
tic
if nargin == 1
    secs_of_data = 1;
end
if nargin < 2
    Ne=800;
end
if nargin < 3
    Ni=200;
end
if nargin < 4
    neuron_params = 'default';
end
if nargin < 5
    fq = []; % induce a background frequency across all neurons.
    % if you want to pass in a frequency, use a tuple: the first number is
    % the fq in Hz the second number is how strong (8 would be very strong).
end
rng(123) % reset the seed. Keep it the same
switch neuron_params
    case 'default'
        %% Excitatory neurons    Inhibitory neurons
        re=rand(Ne,1);       ri=rand(Ni,1);
        a=[0.02*ones(Ne,1);     0.02+0.08*ri];
        b=[0.2*ones(Ne,1);      0.25-0.05*ri]; % induces oscllations if high
        c=[-65+15*re.^2;        -65*ones(Ni,1)];
        d=[8-6*re.^2;           2*ones(Ni,1)];
        thale = 5; thali = 2;
        S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];
    case 'low_firing_rate'
        re=rand(Ne,1);          ri=rand(Ni,1); 
        a=[0.02*ones(Ne,1);     0.02+0.08*ri];
        b=[0.21*ones(Ne,1);      0.25-0.05*ri];
        c=[-65+32*re.^2;        -65*ones(Ni,1)]; % this can change the burstiness
        d=[8-6*re.^2;           2*ones(Ni,1)];
        thale = 3.1; thali = 2.8;
        S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];
    otherwise
        error('wrong parameter type')
end
% induce an oscillation or two.
if isempty (fq)
    OSC = [];
else
    OSC = zeros(1,length(1:secs_of_data*1000));
    for ii = 1:size(fq,1)
        OSC = OSC + sin(linspace(0,2*pi*secs_of_data*fq(ii,1),length(1:secs_of_data*1000 ))) * fq(ii,2);
    end
end
v=-65*ones(Ne+Ni,1);    % Initial values of v
th = 30; % threshold for firing.
%  th = 30 + [randn(Ne,1)*15; randn(Ni,1)*5+1]; % this may not matter.
u=b.*v;                 % Initial values of u
firings=[];             % spike timings
first_ix = 1;

for t=1:secs_of_data*1000            % simulation of 1000 ms
    I = [thale*randn(Ne,1);thali*randn(Ni,1)]; % thalamic input
    if ~isempty(OSC)
        I = I + OSC(t);
    end
    fired=find(v>=th);    % indices of spikes
    if 0
        firings=[firings; t+0*fired,fired];
    else
        if ~isempty(fired)
            end_ix = (first_ix + size(fired,1)-1);
            firings(first_ix:end_ix,:) =  [t+0*fired,fired];
            first_ix = end_ix + 1;
        end
    end
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);
    I=I+sum(S(:,fired),2);
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u);                 % stability
end;
%
if nargout ==0 
    cnt = hist(firings(:,2),1:(Ne + Ni));
    fr = cnt/secs_of_data;
    figure
    subplot(4,1,1:2)
    plot(firings(:,1)/1000,firings(:,2),'.');
    ylabel('Neuron ID')
    subplot(4,1,3)
    r = 1:200:secs_of_data*1000;
    H = histc(firings(:,1),r);
    bar(r/1000,H);axis tight
    ylabel('nSpikes')
    xlabel('s')
    subplot(4,1,4)
    bar(fr);
    axis tight
    ylabel('Hz')
    xlabel('Neuron ID')
    
    nC = round((Ne + Ni)/5);
    nR = round((Ne + Ni)/nC);
    figure
    for ii = 1:(Ne + Ni)
        t = firings(firings(:,2)==ii,1);
        subplot(nR,nC,ii)
        [c,x] = AutoCorr(t*10,5,60);
        plot(x,c,'k'); axis tight;
        box off
    end
end
toc