clear;clc

connfit_tag = 'pls';
config = [];
config.names_MRImetrics_all = {'Ne','Lp','Neloc'}; 
config.which_metrics = [1 2 3];

flag_write_PLSstats= 1;
flag_write_generank = 1;

config.RPKM_thrmin = NaN;        
config.flag_zscoreX = 'zscoreX';  
config.regionlist_touse = [1:11];
config.num_perm_forPLSstats = 10000;   
config.PCA = 0;   
config.PCA_comp = 1; 

config.num_comp_forbs = 1;
config.if_save_comp_vec_forbs = [1]; 
config.if_flip_XS_vec_forbs = [1];   
config.bootnum_forbs = 10000;   

folder_integration = 'Data';

config.agepcd_vec = 3186;
config.list_subjstr = {'pcd3186'};

config.folder_RPKM = 'Data';
config.prefix_RPKMfile = 'RPKM_Model';
config.postfix_RPKMfile = '11NCX';

config.tmetrics = [];
config.folder_GRETNA_myresults = folder_integration;
names_BrainSpan_regions_startsMFC = {'MFC','OFC','DFC','VFC','M1C','S1C','IPC','A1C','STC','ITC','V1C'}; 
for i_metric = 1:length(config.names_MRImetrics_all)
    str_metric = config.names_MRImetrics_all{i_metric};
    tmetrics_full = readtable([config.folder_GRETNA_myresults,filesep,str_metric,'_11ROI_fit_8Years.txt']);
    
    idx_3186 = find(tmetrics_full.agepcd_todo == config.agepcd_vec);
    if length(idx_3186) ~= 1
        error('expected exactly one row with agepcd_todo == 3186 in tmetrics for metric %s', str_metric);
    end
    config.tmetrics.(str_metric) = tmetrics_full(idx_3186, :);
    config.tmetrics.(str_metric).agepcd_todo = [];
    
    if min(cellfun(@isequal,names_BrainSpan_regions_startsMFC, config.tmetrics.(str_metric).Properties.VariableNames)) == 0 
        error('the order of ROIstr in tmetrics is not the same as expected one!');
    end
end

config.folder_integ = [folder_integration,filesep,'PLSstats'];

config.num_group  = length(config.agepcd_vec);
%%
prompt = ['Is this the folder for saving? (enter "y" if yes,"N" if NO) ',config.folder_integ,'\n','xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'];
input_str = input(prompt,'s');
if strcmp(input_str,'y') ~= 1
    error('please check the folder for saving');
end
mkdir(config.folder_integ);

if flag_write_PLSstats == 1
    PLSregression_func(config);
    disp('start to run PLSregression_func.m for PLS stats.....');
end

if flag_write_generank == 1
    PLSregression_BootstrappedWeights_func(config);
    disp('start to run PLSregression_BootstrappedWeights_func.m for ranked gene list.....');
end
