function t =  mat2table(mat, RowNames, VariableNames)
t = table(mat);
t = splitvars(t);
if ~isempty(RowNames)   
    t.Properties.RowNames = RowNames;
end
if nargin==3
    t.Properties.VariableNames = VariableNames;
end
end