%% 2-way ANOVA for SPES across subgroups (e.g. mtle vs. non-mtle)

% close all
% clear all
% clc

load mileage

data = []


nmbpats = 11; % Number of cars from each model, i.e., number of replications
[p,anova_tbl,stats] = anova2(data,nmbpats);

c = multcompare(stats);

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

