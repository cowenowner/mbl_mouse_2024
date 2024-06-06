function [m,max_v] = mode_ksdensity(v)
[f,xi] = ksdensity(v);
[max_v,max_ix] = max(f);
m = xi(max_ix);