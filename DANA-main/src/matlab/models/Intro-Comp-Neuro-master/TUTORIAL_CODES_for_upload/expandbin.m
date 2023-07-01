function [ new_vector ] = expandbin( old_vector, old_dt, new_dt )
%expandbin.m takes a vector of values with very finely spaced time points
%and returns a smaller vector with more coarsely spaced time points
%   [newvector] = expandbin( old_vector, old_dt, new_dt )
%
%   The function requires the following inputs:
%   old_vector contains the values in original time bins.
%   These will be averaged to generate values in the new time bins and
%   returned as new_vector
%
%   old_dt is the original time-step used
%
%   new_dt is the desired new time-step
%
%   When analyzing simulated data this function is useful, because the
%   simulations may require a very small dt for accuracy and stability, yet
%   the results may only need to be analyzed at a less fine resolution.
%
%   This code is needed for Tutorial 3.1 in Chapter 3 of the textbook
%   An Introductory Course in Computational Neuroscience
%   by Paul Miller, Brandeis University, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

old_Nt = length(old_vector);
Nscale = round(new_dt/old_dt);
new_Nt = ceil(old_Nt/Nscale);
new_vector = zeros(1,new_Nt);

for i = 1:new_Nt-1
    new_vector(i) = mean(old_vector((i-1)*Nscale+1:i*Nscale));
end
new_vector(end) = mean(old_vector((new_Nt-1)*Nscale:end));

