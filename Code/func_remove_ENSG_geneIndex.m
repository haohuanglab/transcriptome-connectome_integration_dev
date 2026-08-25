function b = func_remove_ENSG_geneIndex(a)
len = length(a);
idx_split = 16;  
b = cell(len,1);
for i_gene = 1:len
    temp = a{i_gene};
    b{i_gene} = temp(idx_split+1 : end);
end
end