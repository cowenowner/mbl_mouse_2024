function [CI1, CI2] = Vector_Stats_Bootstrap_CImnVM(Ydata, ...
                                B, N, ThetaBarR, KappaHat, alfa)

%Vector_Stats_Bootstrap_CIkappa.m

% Generate confidence interval on concentration parameter of vonMises
%   using resampling (bootstrapping) methods.
% Assume von Mises distribution

% Input variables:
%   Ydata: azimuth data (radians)
%   B: number of bootstrapping iterations
%   N: sample size
%   ThetaBarR: calculated vector mean (radians)
%   KappaHat: calculated kappa estimate
%   alfa: significance level
% Output variables:
%   CI1, CI2:  limits of confidence interval on vector mean

% Called functions:
%   Vector_Stats_Bootstrap_Algor.m
%   VonMises_GenerateRandomRad.m
%   CalcKappa.m

% Ref.: Fisher, 1993, p. 90-91, 199-207, 210-211

% calculate algorithm values for entire sample

[Z0, V0, W0] = Vector_Stats_Bootstrap_Algor(Ydata);

% loop over B iterations

phi = zeros(N,1);
BStraps = zeros(B,1);

for iii = 1:B
    
    % generate random resample
    
    phi = VonMises_GenerateRandomRad(N, ThetaBarR, KappaHat);
%    cc=sum(cos(phi))/N;ss=sum(sin(phi))/N;ss=sqrt(cc*cc+ss*ss)
%    ee=CalcKappa(ss, 0)
     
    % calculate estimate from resample
      
    [Z, V, W] = Vector_Stats_Bootstrap_Algor(phi);
    
    CBSB = Z0 - V0*W*(Z - Z0);
    CB = CBSB(1);
    SB = CBSB(2);
    est = VectMean_arctan(SB, CB); 
    BStraps(iii) = est/57.3;
    
end

% calculate confidence interval

gamma = BStraps - ThetaBarR;
for iii = 1:B
    a = gamma(iii);
    if a < -pi
        gamma(iii) = a + pi
    end
    a = gamma(iii);
    if a > pi
        gamma(iii) = a - pi;
    end
end
    
Sgamma = sort(gamma);
    
ind1 = fix(0.5 + B*alfa/2);
ind2 = B - ind1;
    
CI1 = ThetaBarR + Sgamma(ind1+1);
CI2 = ThetaBarR + Sgamma(ind2);

pi2 = 2*pi;
if CI1 < 0 
    CI1 = CI1 + pi2;
end
if CI1 > pi2
    CI1 = CI1 - pi2;
end
if CI2 < 0 
    CI2 = CI2 + pi2;
end
if CI2 > pi2
    CI2 = CI2 - pi2;
end

clear phi BStraps gamma Sgamma
  
   