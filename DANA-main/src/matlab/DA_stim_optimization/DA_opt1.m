function [objective]  = DA_opt1(STIM, TGT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare similarity to each element in target.
% Dopamine is released depending on how close the pattern is to the
% best match in the target library.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DA_tmp = zeros(size(TGT,1),1);
% Constrain STIM to be positive since it represents the positive ISI from
% the previous stim. Actually, this might not be necessary. A negative ISI
% just means that spike comes before the previous one so it is
% interpretable I think.
% STIM = abs(STIM);
% STIM(STIM<0.0001) = 0;
% TODO: Using a loop is not optimal - could be simple matrix math. For
% example a = STIM-TGT; DA_tmp = -1*a*a';
for jj = 1:size(TGT,1)
    % Pick an objective function to minimize. -1 to force minimization.
    % SSE is the most interpretable right now.
    % 
    %           c = corrcoef(TGT(jj,:),STIM); % r values if more interested in
    %      pattern rather than magnitude.
    %           DA_tmp(jj) = c(2);
    %     DA_tmp(jj)  = dot(TGT(jj,:),STIM); % dot
    %          DA_tmp(jj) = -1*mean(abs(TGT(jj,:)-STIM)); % Abs dist. The
    %          DA_tmp(jj) = 1/mean(abs(TGT(jj,:)-STIM)); % Abs dist. The
    %     smaller the distance, the more dopamine is released.
    %     DA_tmp(jj) = -1*mean((STIM-TGT(jj,:)).^2); % MSE. The smaller the MSE, the more dopamine released.
     DA_tmp(jj) = -1*sum((STIM-TGT(jj,:)).^2); % SSE. faster
% incorporate the firs derivative as well - differences between adjacent
% intervals.
%     DA_tmp(jj) = -1*sum((STIM-TGT(jj,:)).^2) + -.7*sum((diff(STIM)-diff(TGT(jj,:))).^2); % SSE. faster
% Another alternative would be to use a function that optimizes the
% distribution itself - the fit to the histogram - try different firing
% distributions: Pink noise, gaussian, etc... 
% Another idea is to generate the ISIs using a function like a polynomial
% and then just modify the polynomial coeffieicnts instead of each
% individual ISI. That could reduce the search problem and produce more
% generalizable results. Process would be to fit a polynomial to the tgt -
% say of 20 pulses into a 5th order polynomial wiht 5 coefficeints. THen
% fit a polynomial to the stim. Measure the similarity between the stim
% coeffieicnts and the tgt adn that becomes the objective function.

end

% Should add a penalty for negative ISIs.
negpen = double(sum(STIM<0.001))*1.9;

objective = -1*max(DA_tmp)+negpen; % the energy function is maximal dopmaine release. Since optimization needs to minimize, then we have to multiply by -1
% To encourage convergence on one pattern, it could 