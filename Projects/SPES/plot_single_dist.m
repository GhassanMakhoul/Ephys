function f = plot_single_dist(f,A,B,C,distance_for_bar_plot,title_str,x_label,y_label,ymin,ymax,soz_color, pz_color, non_color, jitterAmount)

clf
f.Position = [100 100 540 900];
hold on
soz_bar_data_IN = A(:,distance_for_bar_plot);
pz_bar_data_IN = B(:,distance_for_bar_plot);
non_bar_data_IN = C(:,distance_for_bar_plot);
soz_data_ci95 = 1.96 * std(soz_bar_data_IN, "omitnan")/sqrt(length(soz_bar_data_IN));
pz_data_ci95 = 1.96 * std(pz_bar_data_IN, 'omitnan')/sqrt(length(pz_bar_data_IN));
non_data_ci95 = 1.96 * std(non_bar_data_IN, 'omitnan')/sqrt(length(non_bar_data_IN));
b=bar([mean(soz_bar_data_IN, 'omitnan'), mean(pz_bar_data_IN,'omitnan'), mean(non_bar_data_IN, 'omitnan')]);
b.FaceColor = 'flat';
b.CData(1,:) = soz_color; % SOZ
b.CData(2,:) = pz_color; % PZ
b.CData(3,:) = non_color; % Non
%Dummy objects for legend use
bh(1) = bar(nan,nan,'FaceColor',soz_color);
bh(2) = bar(nan,nan,'FaceColor',pz_color);
bh(3) = bar(nan,nan,'FaceColor',non_color);
legend(bh,'SOZ','PZ','Non','AutoUpdate','off')
errorbar([mean(soz_bar_data_IN, 'omitnan'), mean(pz_bar_data_IN ,'omitnan'), mean(non_bar_data_IN,'omitnan')],...
    [soz_data_ci95, pz_data_ci95, non_data_ci95],'LineStyle','none','Color','k','LineWidth',2)
scatter(ones(length(soz_bar_data_IN),1), soz_bar_data_IN,60,'MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerEdgeColor','k','LineWidth',1,'jitter', 'on', 'jitterAmount', jitterAmount)
scatter(ones(length(soz_bar_data_IN),1)*2, pz_bar_data_IN,60,'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','k','LineWidth',1,'jitter', 'on', 'jitterAmount', jitterAmount)
scatter(ones(length(soz_bar_data_IN),1)*3, non_bar_data_IN,60,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor','k','LineWidth',1,'jitter', 'on', 'jitterAmount', jitterAmount)
title(title_str);
ylim([ymin ymax])
ylabel(y_label)
hold off
% Stats on IN bar chart
[P_IN,ANOVATAB_IN,STATS_IN] = anova1([soz_bar_data_IN, pz_bar_data_IN, non_bar_data_IN],[],'off');
ANOVATAB_IN
annotation('textbox', [0.15, 0.058, 0.1, 0.1], 'String', sprintf("One-way ANOVA p=%0.4f",P_IN))
[c,m,h,nms] = multcompare(STATS_IN,'display','off');
annotation('textbox', [0.62, 0.09, 0.1, 0.1], 'String', sprintf("SOZ-PZ p=%0.4f\nSOZ-Non p=%0.4f\nPZ-Non p=%0.4f",c(1,end),c(2,end),c(3,end)))

end