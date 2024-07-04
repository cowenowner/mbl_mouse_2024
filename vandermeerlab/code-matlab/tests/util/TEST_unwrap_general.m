%% setup

%% basic test
in = [9 10 1];
range = [0 10];

out = unwrap_general(in, range(1), range(2));
assert(isequal(out, [0 1 2]));
