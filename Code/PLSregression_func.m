function PLSregression_func(config)
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
num_perm                    = config.num_perm_forPLSstats;   
if_pca                      = config.PCA;  
num_components              = config.PCA_comp; 
if ~isnumeric(RPKM_thrmin); error('RPKM_thrmin is not numeric! If not want to thr, enter 0'); end

%% 
gpuDevice();
names_BrainSpan_regions_startsMFC = {'MFC','OFC','DFC','VFC','M1C','S1C','IPC','A1C','STC','ITC','V1C'}; 

load('Data\gene_namelist_60155genes.mat');  
gene_namelist_all = gene_namelist;

num_comp_full = length(names_BrainSpan_regions_startsMFC) - 1;  

if_exclude_RPKMabove1000 = 0;  

fn_DEXlist = ['Data',filesep,'GenesSpatialDEX_geneidxlistOutof60155_protein_coding.mat']; 
DEXlist = load(fn_DEXlist); DEXlist = DEXlist.geneidxlist_DEX_sorted;

disp(['================================================= which_metrics = ',num2str(which_metrics),'=================================================']);
    Tresults = []; 
    Tresults.p_pls = zeros(num_group, num_comp_full);  
    Tresults.p_boot = zeros(num_group, num_comp_full);  
    Tresults.ExplVarinY = zeros(num_group, num_comp_full);  
    Tresults.p_corr = cell(num_group, num_comp_full);  
    Tresults.r_corr = cell(num_group, num_comp_full);  
    Tresults.p_perm_corr = cell(num_group, num_comp_full);  
    Tresults.XS = cell(num_group, num_comp_full); 
    Tresults.geneidxlistOutof60155 =  cell(num_group, 1); 
tic;

i_group = 1;
    disp(['--- ','which group = ',num2str(i_group),' ---']);
    if size(which_metrics,1)~=1; error('enter a row vector of which_metrics!'); end 

    if_debugPLS = 0;  

    if length(regionlist_touse) < length(names_BrainSpan_regions_startsMFC)
        warning('Not all regions will be used for PLS regression!');
        list_missing_region = setdiff(1:length(names_BrainSpan_regions_startsMFC), regionlist_touse);
        for i_region = list_missing_region
            disp(['missing: ', names_BrainSpan_regions_startsMFC{i_region}]);
        end
    end

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

    if size(X,2)~= 60155
        error('should not use this geneidx list to extract DEX genes!');
    end

    if isnan(RPKM_thrmin)   
        geneidxlist_minRPKM1inNCX = DEXlist;   
    else
        geneidxlist_minRPKM1inNCX = intersect(find(min(X,[],1)>RPKM_thrmin)',...  
                                            DEXlist);                            
    end
    disp(['#=', num2str(length(geneidxlist_minRPKM1inNCX)),' genes will be used...']);
    num_metric = size(Y, 2);
    X = X(:,geneidxlist_minRPKM1inNCX);

    X = X(regionlist_touse,:);  

    Y = Y(regionlist_touse,:);   

    if strcmp(flag_zscoreX,'zscoreX')
        X = zscore(X);   
    elseif strcmp(flag_zscoreX, 'NozscoreX')
        X = X;
    else
        error('invalid flag_zscoreX!');
    end

    Y_zscore = zscore(Y);  
    num_metric = size(Y_zscore,2);
    if if_pca
        if num_components == 0
            num_components = min(size(score, 2), num_comp_full);
        end
        Y_pca = zeros(num_NCX, num_components);
        [coeff, score, ~, ~, explained, ~] = pca(Y_zscore);
        Y_pca = score(:, 1:num_components);
        Y_zscore = Y_pca;
    end

    if if_debugPLS
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,statsPLS] = plsregress_CZtry(X, Y_zscore, num_comp_full);  
    else
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,statsPLS] = plsregress(X, Y_zscore, num_comp_full, "options",statset("UseParallel",true)); 
    end
    wPLS = statsPLS.W;

    [r_CompxMetric, p_CompxMetric] = corr(XS, Y_zscore);
    Rsq_CompxMetric = r_CompxMetric.* r_CompxMetric;
    PCTVARinY_found_vec = cumsum(PCTVAR(2,:));    % in Y
    
    num_boot = num_perm;
    which_perm = 'X';   
    num_metric = size(Y_zscore,2);
    perm_PCTVARinY_vec = zeros(num_perm, num_comp_full);
    perm_Rsq_CompxMetric = zeros(num_comp_full, num_metric, num_perm); 
    boot_PCTVARinY_vec = zeros(num_boot, num_comp_full);   

    for i_perm = 1:num_perm
        if strcmp(which_perm, 'X')
            order_xperm = randperm(size(X, 1));
            X_perm = X(order_xperm,:);
        else
            order_perm = randperm(size(Y_zscore,1));
            Y_perm = Y_zscore(order_perm,:);
        end
        num_metric = size(Y_zscore,2);
        
        if strcmp(which_perm, 'X')
            [~,~,XS_perm,YS_perm,MSE_perm,PCTVAR_perm] = plsregress(X_perm, Y_zscore, num_comp_full, "options",statset("UseParallel",true));
        else
            [~,~,XS_perm,YS_perm,MSE_perm,PCTVAR_perm] = plsregress(X, Y_perm, num_comp_full, "options",statset("UseParallel",true));
        end
        
        perm_PCTVARinY_vec(i_perm, :) = cumsum(PCTVAR_perm(2,:));  
        if strcmp(which_perm, 'X')
            perm_Rsq_CompxMetric(:,:,i_perm) = (corr(XS_perm, Y_zscore)).^2;
        else
            perm_Rsq_CompxMetric(:,:,i_perm) = (corr(XS_perm, Y_perm)).^2;
        end  
    end

    p_explvar_vec = zeros(1, num_comp_full);
    p_boot_vec = zeros(1, num_comp_full);
    p_perm_corr = zeros(num_metric, num_comp_full);   % metric x ncomp

    for i_comp = 1:num_comp_full
        PCTVARinY_found = PCTVARinY_found_vec(i_comp);
        perm_PCTVARinY = perm_PCTVARinY_vec(:,i_comp);
        boot_PCTVARinY = boot_PCTVARinY_vec(:,i_comp);
        p_explvar_vec(1,i_comp) = length(find(perm_PCTVARinY>=PCTVARinY_found))./num_perm;
        p_boot_vec(1,i_comp) = length(find(boot_PCTVARinY>=PCTVARinY_found))./num_boot;
        for i_metric = 1:num_metric
            Rsq_corr_found = Rsq_CompxMetric(i_comp, i_metric);
            perm_Rsq_corr_vec = squeeze(perm_Rsq_CompxMetric(i_comp, i_metric, :));
            p_perm_corr(i_metric, i_comp) = length(find(perm_Rsq_corr_vec >= Rsq_corr_found )) ./ num_perm;
        end
    end

    Tresults.p_pls(i_group, :) = p_explvar_vec;
    Tresults.p_boot(i_group, :) = p_boot_vec;
    Tresults.ExplVarinY(i_group, :) = PCTVARinY_found_vec;
    for i_comp = 1:num_comp_full
        Tresults.p_perm_corr{i_group, i_comp} = p_perm_corr(:, i_comp);
        Tresults.r_corr{i_group, i_comp} = r_CompxMetric(i_comp, :)';
        Tresults.p_corr{i_group, i_comp} = p_CompxMetric(i_comp, :)';
        Tresults.XS{i_group, i_comp} = XS(:,i_comp);
    end
    Tresults.geneidxlistOutof60155{i_group,1} = geneidxlist_minRPKM1inNCX;  

    Tresults.p_pls = mat2table(Tresults.p_pls, list_subjstr); 
    Tresults.p_boot = mat2table(Tresults.p_boot, list_subjstr);
    Tresults.ExplVarinY = mat2table(Tresults.ExplVarinY, list_subjstr);
    Tresults.p_corr = mat2table(Tresults.p_corr,list_subjstr);
    Tresults.r_corr = mat2table(Tresults.r_corr, list_subjstr);
    Tresults.p_perm_corr = mat2table(Tresults.p_perm_corr, list_subjstr);
    Tresults.XS         = mat2table(Tresults.XS,         list_subjstr);
    Tresults.geneidxlistOutof60155 = mat2table(Tresults.geneidxlistOutof60155, list_subjstr);

    Tresults.ROIstr = names_BrainSpan_regions_startsMFC(regionlist_touse);
    Tresults.config = config;

    temp = dbstack;
    Tresults.functionName = temp.name;

    temp = which_metrics_str_generator(which_metrics, names_MRImetrics_all); temp = temp(2:end);
    if if_pca
        fn_Tresults = [folder_integ, filesep, temp, '_PCA', num2str(num_components),'.mat'];
    else
        fn_Tresults = [folder_integ, filesep, temp,'.mat'];
    end
    save(fn_Tresults,'Tresults');
toc; 

end   