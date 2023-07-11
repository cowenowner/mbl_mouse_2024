%
final_seq = [];
for i = 1:350
    toss = rand();
    if toss < 0.9
        final_seq = [final_seq, mod(i,3)+1];
    else
        temp = randperm(3);
        final_seq = [final_seq, temp(1)];
    end
end