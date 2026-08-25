function func_ReadGeneid_SaveasTxt(fn)
T = readtable(fn);
PLS = T.Geneid;
PLS_geneName = func_remove_ENSG_geneIndex(PLS);

[~,~,ext] = fileparts(fn);
if isempty(ext)
    error('there is no extension in filename!! Please add .csv at the end of filename!');
end
fn_out = strrep(fn,ext,'.txt');  
fid = fopen(fn_out,'w');
fprintf(fid, '%s\n',PLS_geneName{:});
fclose(fid);
end