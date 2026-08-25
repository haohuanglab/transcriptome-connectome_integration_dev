function [avg_diff, all_diff] = func_diff_rank_in_list(list1,list2)
len = length(list1);
if length(list1) ~= length(list2)
    error('length of two lists are not the same!');
end

all_diff = zeros(len,1);
for i = 1:len
    a = list1(i);
    i_inlist2 = find(list2==a);
    all_diff(i) = abs(i - i_inlist2) ./ len;
end
avg_diff = mean(all_diff);
end