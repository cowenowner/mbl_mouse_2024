% Montague model.
% Montague PR, McClure SM, Baldwin PR, Phillips PEM, Budygin E a, Stuber GD, et al. Dynamic gain control of dopamine delivery in freely moving animals. The Journal of Neuroscience : The Official Journal of the Society for Neuroscience. 2004;24:1754–1759.
% % NOTE _ I also wrote a python version in the python folder. 
%% Use the values from a given row in Table 1.
%%
% 3 exp function/dynamical equation model: short facilitatio, short depression,
% long depression.
clearvars
IBI_s = 1;
a_0_nM = 10;
tau_s(1) = 4.16;
k(1) = 1.012; % Facilitation because >1
tau_s(2) = 3.53;
k(2) = .997; % Depression because <1
dt_s = 0.01;
Vm_uM_sec = 4;
Km_uM = .2;
I_init = [1 1]; % If it's 1, then no slow change - seems pointless then.
A_init = .5;
C_init = 0.001;
%
n_factors = length(tau_s);
tmax_s = 20;
t_s = 0:dt_s:tmax_s;
spikes = rand(size(t_s)) > 0.95;
% Ij is a hiddent dymaic process. we modeled the concentration of dopamine
% added by each impulse as a product of two time-dependent functions
% p(t)A(t). A(t) is a function that depends on independent facilitation and
% depression components Ij(t), and its initial value, a0 . p(t) is a
% function that models the exact impulse times for a specific pattern of
% impulses. As detailed below, each Ij(t) possesses "kick and relax"
% dynamics.
I = zeros(length(t_s),n_factors);
I(1,:) = I_init;
C = zeros(size(t_s));
A = zeros(size(t_s));
C(1) = C_init;
A(1) = A_init;
V = zeros(1,n_factors);

for ii = 2:length(t_s)
    % dCdt = r*Cp - Vm/(1+Km/C);
    % Compute teh 'kick' factor * k_1
    % spikes(ii) = s = 0 = no spike, s = 1 = spike.
    for jj = 1:n_factors
        I(ii,jj) = I(ii-1,jj) + dt_s*(1-I(ii-1,jj))/tau_s(jj); % here is where I get confused - should it be I - 1 or 1- I - something gets flipped in the conversion the euler method I think
    end

    if spikes(ii)
        for jj = 1:length(V)
            V(jj) = k(jj)*I(ii,jj);
        end
        A(ii) = a_0_nM*prod(V);
        %         disp('spke')
    end

    C(ii) = C(ii-1) + dt_s*(spikes(ii)*A(ii) - Vm_uM_sec/(1 + Km_uM/C(ii-1)));

    %
end

figure
subplot(2,1,1)
plot(t_s, I(:,1),t_s, I(:,2))
ylabel('I')
subplot(2,1,2)
plot(t_s, C)
hold on
plot(t_s(spikes),t_s(spikes)*0 ,'r+')
axis tight
xlabel('sec');ylabel('C')







