function [p_CompxMetric, slope_CompxMetric] = func_corr_metric_PLSscores(XS, Y)
num_comp = size(XS,2);
num_metric = size(Y,2);

p_CompxMetric = zeros(num_comp,num_metric);
Rsq_CompxMetric = zeros(num_comp,num_metric);
slope_CompxMetric = zeros(num_comp,num_metric);
intercept_CompxMetric = zeros(num_comp,num_metric);
for i_metric = 1:num_metric
    for i_comp = 1:num_comp
        x = XS(:,i_comp);
        y = Y(:,i_metric);
        [b,~,~,~,stats] = regress(y,[ones(length(x),1),x]);
        p_CompxMetric(i_comp,i_metric) = stats(3);
        slope_CompxMetric(i_comp,i_metric) = b(2);
        intercept_CompxMetric(i_comp,i_metric) = b(1);
        Rsq_CompxMetric(i_comp,i_metric) = stats(1);
    end
end