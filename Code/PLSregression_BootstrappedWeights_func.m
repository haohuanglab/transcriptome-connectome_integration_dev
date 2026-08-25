function [PCTVARinY_found_vec] = PLSregression_BootstrappedWeights_func(config)
agepcd_vec                  = config.agepcd_vec;
list_subjstr                = config.list_subjstr;   
tmetrics                    = config.tmetrics;
folder_integ                = config.folder_integ;
names_MRImetrics_all        = config.names_MRImetrics_all;
which_metrics               = config.which_metrics;  
num_group                   = config.num_group;
folder_RPKM                 = config.folder_RPKM;
prefix_RPKMfile             = config.prefix_RPKMfile;
postfix_RPKMfile            = config.postfix_RPKMfile;
RPKM_thrmin                 = config.RPKM_thrmin;   
flag_zscoreX                = config.flag_zscoreX;    
regionlist_touse            = config.regionlist_touse;   

num_comp                    = config.num_comp_forbs;  
bootnum                     = config.bootnum_forbs;
if_save_comp_vec            = config.if_save_comp_vec_forbs;    
if_flip_XS_vec              = config.if_flip_XS_vec_forbs;  

if ~isnumeric(RPKM_thrmin); error('RPKM_thrmin is not numeric! If not want to thr, enter 0'); end
if num_comp > 3
    error('current PLSregression_BootstrappedWeights_func.m does not accept num_comp > 3!');
end
if length(if_flip_XS_vec) ~= num_comp
    error('length of vector "if_flip_XS_vec" should be equal to "num_comp_forbs" ');
end
if length(if_save_comp_vec) ~= num_comp
    error('length of vector "if_save_comp_vec" should be equal to "num_comp_forbs" ');
end
%% 
names_BrainSpan_regions_startsMFC = {'MFC','OFC','DFC','VFC','M1C','S1C','IPC','A1C','STC','ITC','V1C'}; 

gene_namelist_all = load('Data\gene_namelist_60155genes.mat'); 
gene_namelist_all = gene_namelist_all.gene_namelist;

fn_DEXlist = ['Data',filesep,'GenesSpatialDEX_geneidxlistOutof60155_protein_coding.mat'];
DEXlist = load(fn_DEXlist); DEXlist = DEXlist.geneidxlist_DEX_sorted;

if length(regionlist_touse) < length(names_BrainSpan_regions_startsMFC)
    warning('Not all regions will be used for PLS regression!');
    list_missing_region = setdiff(1:length(names_BrainSpan_regions_startsMFC), regionlist_touse);
    for i_region = list_missing_region
        disp(['missing: ', names_BrainSpan_regions_startsMFC{i_region}]);
    end
end
%% 
disp(['=========================================  which_metrics = ',num2str(which_metrics),'=================================================']);
    mylog = [];
    mylog.ExplVarinY = zeros(num_group, num_comp);  
    mylog.num_flips_XS_bootstrap = zeros(num_group, num_comp);
    mylog.rankedgenelist_afterbs_avgdiff = zeros(num_group, num_comp);
    mylog.rankedgenelist_afterbs_perc_diffnolargerthan0p05 = zeros(num_group, num_comp);  
    mylog.XS_orig = cell(num_group, num_comp); 
    mylog.XS_used = cell(num_group, num_comp); 
    mylog.geneidxlistOutof60155 = cell(num_group,1);  
    mylog.if_flipweight_bs_by_numcomp = zeros(bootnum, num_comp);
tic;

i_group = 1;
disp(['--- ','which group = ',num2str(i_group),' ---']);
if size(which_metrics,1)~=1; error('enter a row vector of which_metrics!'); end 
 
num_total_genes = 60155;
num_NCX = 11;

agepcd = 3186;   
fn_RPKM = [folder_RPKM, filesep,prefix_RPKMfile, '_pcd',num2str(agepcd),'_',postfix_RPKMfile,'.txt'];
genemat = [];
genemat = readtable(fn_RPKM);
if ~isequal(names_BrainSpan_regions_startsMFC, genemat.Properties.VariableNames)
    error('The order of region names from RPKM.txt is not aligned to names in Li 2018');
end
if height(genemat)~= num_total_genes
    warning(['number of genes in genemat is not ',num2str(num_total_genes)]);
end

X = genemat{:,:}';

 
Y = zeros(num_NCX,length(which_metrics));
for i_metric = 1:length(which_metrics)
    idx_metric = which_metrics(i_metric);
    str_metric = names_MRImetrics_all{idx_metric};
    Y(:,i_metric) = (tmetrics.(str_metric){i_group,:})';
end

%% 
if size(X,2)~= 60155
    error('should not use this geneidx list to extract DEX genes!');
end

if isnan(RPKM_thrmin)  
    geneindex = DEXlist;   
else
    geneindex = intersect(find(min(X,[],1)>RPKM_thrmin)',... 
                          DEXlist);                            
end
disp(['#=', num2str(length(geneindex)),' genes will be used...']);


X = X(:,geneindex);
gene_namelist = gene_namelist_all(geneindex); 

X = X(regionlist_touse,:);  
Y = Y(regionlist_touse,:); 

if strcmp(flag_zscoreX,'zscoreX')
    X = zscore(X);   
elseif strcmp(flag_zscoreX, 'NozscoreX')
    X = X;
else
    error('invalid flag_zscoreX!');
end

Y_orig = Y;
Y = zscore(Y);   

[XL,YL,XS,YS,BETA,PCTVAR,~,statsPLS]=plsregress(X,Y,num_comp);
PCTVARinY_found_vec = cumsum(PCTVAR(2,1:num_comp));
R1=corr(XS(:,1),Y_orig);   
if num_comp>1
    R2 = corr(XS(:,2),Y_orig);
end
if ~isempty(find(corr(XS,Y)<0)); warning('needs to flip XS!'); end     

[p_CompxMetric, slope_CompxMetric] = func_corr_metric_PLSscores(XS, Y);  

XS_orig = XS;  

if if_flip_XS_vec(1) == 1   
    disp('flipping XS 1 !');
    statsPLS.W(:,1)=-1*statsPLS.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if num_comp>=2 && if_flip_XS_vec(2) == 1  
    disp('flipping XS 2 !');
    statsPLS.W(:,2)=-1*statsPLS.W(:,2);
    XS(:,2)=-1*XS(:,2);
end
if num_comp>=3 && if_flip_XS_vec(3) == 1  
    disp('flipping XS 3 !');
    statsPLS.W(:,3) = -1*statsPLS.W(:,3);
    XS(:,3) = -1*XS(:,3);
end
[PLS1w,x1] = sort(statsPLS.W(:,1),'descend');
PLS1ids=gene_namelist(x1);  
geneindex1=geneindex(x1);
if num_comp>=2
    [PLS2w,x2] = sort(statsPLS.W(:,2),'descend');
    PLS2ids=gene_namelist(x2);
    geneindex2=geneindex(x2);
end
if num_comp>=3
    [PLS3w,x3] = sort(statsPLS.W(:,3),'descend');
    PLS3ids=gene_namelist(x3); 
    geneindex3=geneindex(x3);
end

disp('  Bootstrapping - could take a while')
PLS1weights = zeros(size(X,2),bootnum);
if_flip_vec_XS1 = zeros(bootnum,1);
if num_comp >=2
    PLS2weights = zeros(size(X,2),bootnum);
    if_flip_vec_XS2 = zeros(bootnum,1);
end
if num_comp >=3
    PLS3weights = zeros(size(X,2),bootnum);
    if_flip_vec_XS3 = zeros(bootnum,1);
end

myresample_mat = zeros(size(X,1), bootnum);
for i = 1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);  
    myresample_mat(:, i) = myresample;
    Xr=X(myresample,:); 
    Yr=Y(myresample,:); 
    [~,~,~,~,~,~,~,statsBS]=plsregress(Xr,Yr,num_comp); 
    
    temp=statsBS.W(:,1);
    newW=temp(x1); 
    if corr(PLS1w,newW)<0
        newW=-1*newW;
        if_flip_vec_XS1(i) = 1; 
    end
    PLS1weights(:,i) = newW;
    
    if num_comp >= 2
        temp=statsBS.W(:,2);
        newW=temp(x2);
        if corr(PLS2w,newW)<0 
            newW=-1*newW;
            if_flip_vec_XS2(i) = 1;  
        end
        PLS2weights(:,i)= newW; 
    end
    
    if num_comp >=3
        temp=statsBS.W(:,3);
        newW=temp(x3);
        if corr(PLS3w,newW)<0 
            newW=-1*newW;
            if_flip_vec_XS3(i) = 1;  
        end
        PLS3weights(:,i)= newW;
    end
end

num_flips_XS_bootstrap_vec = zeros(1,num_comp);
num_flips_XS_bootstrap_vec(1) = length(find(if_flip_vec_XS1==1));               
disp(['number of flips in XS1_bootstrap: ', num2str(num_flips_XS_bootstrap_vec(1))])
mylog.if_flipweight_bs_by_numcomp(:,1) = if_flip_vec_XS1;
if num_comp>=2
    num_flips_XS_bootstrap_vec(2) = length(find(if_flip_vec_XS2==1)); 
    disp(['number of flips in XS2_bootstrap: ', num2str(num_flips_XS_bootstrap_vec(2))]); 
    mylog.if_flipweight_bs_by_numcomp(:,2) = if_flip_vec_XS2;
end
if num_comp>=3
    num_flips_XS_bootstrap_vec(3) = length(find(if_flip_vec_XS3==1)); 
    disp(['number of flips in XS3_bootstrap: ', num2str(num_flips_XS_bootstrap_vec(3))]); 
    mylog.if_flipweight_bs_by_numcomp(:,3) = if_flip_vec_XS3;
end

PLS1sw=std(PLS1weights');   
temp1=PLS1w./PLS1sw';
if num_comp>=2
    PLS2sw=std(PLS2weights');
    disp(['PLS2: number of genes with 0 std = ',num2str(length(find(PLS2sw==0)))]);
    temp2=PLS2w./PLS2sw';
end
if num_comp>=3
    PLS3sw=std(PLS3weights');
    disp(['PLS3: number of genes with 0 std = ',num2str(length(find(PLS3sw==0)))]);
    temp3=PLS3w./PLS3sw';
end

[Z1,ind1]=sort(temp1,'descend');  
PLS1=PLS1ids(ind1); 
geneindex1_afterBS=geneindex1(ind1); 
if num_comp>=2
    [Z2,ind2]=sort(temp2,'descend');
    PLS2=PLS2ids(ind2);
    geneindex2_afterBS=geneindex2(ind2);
end
if num_comp>=3
    [Z3,ind3]=sort(temp3,'descend');
    PLS3=PLS3ids(ind3);
    geneindex3_afterBS=geneindex3(ind3);
end

PLS1_geneName = func_remove_ENSG_geneIndex(PLS1);
if num_comp>=2; PLS2_geneName = func_remove_ENSG_geneIndex(PLS2); end
if num_comp>=3; PLS3_geneName = func_remove_ENSG_geneIndex(PLS3); end

[avg_diff, all_diff] = func_diff_rank_in_list(geneindex1, geneindex1_afterBS); 
mylog.rankedgenelist_afterbs_avgdiff(i_group, 1) = avg_diff;
mylog.rankedgenelist_afterbs_perc_diffnolargerthan0p05(i_group,1) = length(find(all_diff<=0.05))./length(geneindex1);
disp(['PLS1: mean difference=',num2str(avg_diff),', percent of genes with diff percent below 5% = ', num2str(length(find(all_diff<=0.05))./length(geneindex1))]);
if num_comp >=2   
    [avg_diff, all_diff] = func_diff_rank_in_list(geneindex2, geneindex2_afterBS); 
    mylog.rankedgenelist_afterbs_avgdiff(i_group, 2) = avg_diff;
    mylog.rankedgenelist_afterbs_perc_diffnolargerthan0p05(i_group,2) = length(find(all_diff<=0.05))./length(geneindex2);
    disp(['PLS2: mean difference=',num2str(avg_diff),', percent of genes with diff percent below 5% = ', num2str(length(find(all_diff<=0.05))./length(geneindex2))]);
end
if num_comp >=3  
    [avg_diff, all_diff] = func_diff_rank_in_list(geneindex3, geneindex3_afterBS); 
    mylog.rankedgenelist_afterbs_avgdiff(i_group, 3) = avg_diff;
    mylog.rankedgenelist_afterbs_perc_diffnolargerthan0p05(i_group,3) = length(find(all_diff<=0.05))./length(geneindex3);
    disp(['PLS3: mean difference=',num2str(avg_diff),', percent of genes with diff percent below 5% = ', num2str(length(find(all_diff<=0.05))./length(geneindex3))]);
end
 
temp = which_metrics_str_generator(which_metrics, names_MRImetrics_all); 
foldername_tosave = temp(2:end); 
folder_tosave = [folder_integ, filesep,foldername_tosave];
mkdir(folder_tosave);

tosave_str = ['agepcd',num2str(agepcd)];

if if_save_comp_vec(1) == 1   
    fn_out1 = [folder_tosave, filesep,...
              'geneWeightsBS','_',tosave_str,'_PLS1','.csv'];
    fid1 = fopen(fn_out1,'w');
    fprintf(fid1, '%s,%s,%s\n', 'Geneid', 'gene index in all 60155 genes', 'Bootstrapped Weight');
    for i=1:length(gene_namelist)
      fprintf(fid1,'%s,%d,%f\n', PLS1{i}, geneindex1_afterBS(i), Z1(i));
    end
    fclose(fid1);
    
    fn_out1_rev = strrep(fn_out1, '.csv','_rev.csv');
    fid1 = fopen(fn_out1_rev,'w');
    fprintf(fid1, '%s,%s,%s\n', 'Geneid', 'gene index in all 60155 genes', 'Bootstrapped Weight');
    for i=length(gene_namelist): (-1) :1  
      fprintf(fid1,'%s,%d,%f\n', PLS1{i}, geneindex1_afterBS(i), Z1(i));
    end
    fclose(fid1);
    
    func_ReadGeneid_SaveasTxt(fn_out1);  
    func_ReadGeneid_SaveasTxt(fn_out1_rev);  
end

if num_comp>=2 && if_save_comp_vec(2) == 1  
    fn_out2 = [folder_tosave, filesep,...
              'geneWeightsBS_',tosave_str,'_PLS2','.csv'];
    fid2 = fopen(fn_out2,'w');
    fprintf(fid2, '%s,%s,%s\n', 'Geneid', 'gene index in all 60155 genes', 'Bootstrapped Weight');
    for i=1:length(gene_namelist)
      fprintf(fid2,'%s,%d,%f\n', PLS2{i}, geneindex2_afterBS(i), Z2(i));
    end
    fclose(fid2);
    
    fn_out2_rev = strrep(fn_out2, '.csv','_rev.csv');
    fid2 = fopen(fn_out2_rev,'w');
    fprintf(fid2, '%s,%s,%s\n', 'Geneid', 'gene index in all 60155 genes', 'Bootstrapped Weight');
    for i=length(gene_namelist):(-1):1
      fprintf(fid2,'%s,%d,%f\n', PLS2{i}, geneindex2_afterBS(i), Z2(i));
    end
    fclose(fid2);
    
    func_ReadGeneid_SaveasTxt(fn_out2);  
    func_ReadGeneid_SaveasTxt(fn_out2_rev);  
end

if num_comp>=3 && if_save_comp_vec(3) == 1 
    fn_out3 = [folder_tosave, filesep,...
              'geneWeightsBS_',tosave_str,'_PLS3','.csv'];
    fid3 = fopen(fn_out3,'w');
    fprintf(fid3, '%s,%s,%s\n', 'Geneid', 'gene index in all 60155 genes', 'Bootstrapped Weight');
    for i=1:length(gene_namelist)
      fprintf(fid3,'%s,%d,%f\n', PLS3{i}, geneindex3_afterBS(i), Z3(i));
    end
    fclose(fid3);
    
    fn_out3_rev = strrep(fn_out3, '.csv','_rev.csv');
    fid3 = fopen(fn_out3_rev,'w');
    fprintf(fid3, '%s,%s,%s\n', 'Geneid', 'gene index in all 60155 genes', 'Bootstrapped Weight');
    for i=length(gene_namelist):(-1):1
      fprintf(fid3,'%s,%d,%f\n', PLS3{i}, geneindex3_afterBS(i), Z3(i));
    end
    fclose(fid3);
    
    func_ReadGeneid_SaveasTxt(fn_out3);  
    func_ReadGeneid_SaveasTxt(fn_out3_rev);  
end

fn_out = [folder_tosave, filesep,...
          'myresample_matOfRegionIndex_',tosave_str,'.mat'];
save(fn_out, 'myresample_mat');

fclose all;

mylog.ExplVarinY(i_group, :) = PCTVARinY_found_vec;
mylog.num_flips_XS_bootstrap(i_group,:) = num_flips_XS_bootstrap_vec;
for i_comp = 1:num_comp
    mylog.XS_orig{i_group, i_comp} = XS_orig(:,i_comp);
    mylog.XS_used{i_group, i_comp} = XS(:,i_comp);
end
mylog.geneidxlistOutof60155{i_group,1} = geneindex; 

mylog.ExplVarinY                = mat2table(mylog.ExplVarinY, list_subjstr);
mylog.XS_orig                   = mat2table(mylog.XS_orig, list_subjstr);
mylog.XS_used                   = mat2table(mylog.XS_used, list_subjstr);
mylog.geneidxlistOutof60155     = mat2table(mylog.geneidxlistOutof60155, list_subjstr);
mylog.num_flips_XS_bootstrap    = mat2table(mylog.num_flips_XS_bootstrap, list_subjstr);
mylog.rankedgenelist_afterbs_avgdiff    = mat2table(mylog.rankedgenelist_afterbs_avgdiff, list_subjstr);
mylog.rankedgenelist_afterbs_perc_diffnolargerthan0p05  = mat2table(mylog.rankedgenelist_afterbs_perc_diffnolargerthan0p05, list_subjstr);

mylog.ROIstr = names_BrainSpan_regions_startsMFC(regionlist_touse);
mylog.config = config;

temp = dbstack;
mylog.functionName = temp.name;

fn_mylog = [folder_tosave, filesep, 'mylog.mat'];
save(fn_mylog, 'mylog');

toc;
end  
   


