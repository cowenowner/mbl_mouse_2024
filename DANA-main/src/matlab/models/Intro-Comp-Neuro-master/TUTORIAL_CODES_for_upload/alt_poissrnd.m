function [ num_events ] = alt_poissrnd( rate_dt )
%alt_poissrnd Produces random integers from a Poisson distribution with a
%parameter, lambda = rate_dt that can vary from time bin to time bin.
%
%   rate_dt is an input vector of a time-varying or constant rate,
%   multiplied by the time-step, dt, such that the corresponding bin has on
%   average rate_dt spikes. 
%
%   Output is a vector of integers of the same size as the input vector. 
%   Each element of the output is chosen from the Poisson distribution with
%   paramter given by that element of the input vector.
%
%   This code is required for Tutorial 3.2 (Part B) of Chapter 3 of the
%   textbook,
%   An Introductory Course in Computational Neuroscience
%   by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    num_events = zeros(size(rate_dt));  % initialize with no events per bin
    exprdt = exp(-rate_dt);             % calculated once and used frequently
    rnums = rand(size(rate_dt));        % random numbers to test probabilities
%
    sum_probs = exprdt;                 % cumulative sum of all probabilities
    this_prob = exprdt;                 % probability of N events
    indices = 1:length(rate_dt);                     % bins to test
    N = 0;                              % number of events, N, in a bin
 %
    while ( ~isempty(indices) )         % while there are bins still to test
       indices = find(rnums > sum_probs);   % these bins have more than N events
       N = N+ 1;               % increase N by 1.
       num_events(indices) = N;     % add the new N to each bin
       this_prob = this_prob.*(rate_dt)/N;    % Poisson prob of new N
       sum_probs = sum_probs + this_prob;           % cumulative sum of probs
    end

end

