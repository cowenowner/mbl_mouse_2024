function [O, start_end_times] = Bin_IMU_data(data, edges)
%%%%%%%%
% INPUT: data needs to be a two coloumn matrix with first col of timestamps and second col as the data; edges is a vector of bin edge times 

time = data(:,1);
% Convert the edges into start and end times
tmp1 = edges(1:end-1);
tmp2 = edges(2:end)-eps;
start_end_times = [tmp1(:) tmp2(:)];

% Bin the IMU data 
O = [];
for jj =1:length(start_end_times)
    x = time(:,1) > start_end_times(jj,1) & time(:,1) < start_end_times(jj,2);
    O(jj,1) = mean(data(x,2));
end