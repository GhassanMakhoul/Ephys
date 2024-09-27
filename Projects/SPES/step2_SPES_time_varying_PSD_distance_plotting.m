%% SPES Distance PSD

% Use the generated SPES structs for each patient to plot specrl plots for
% varying distance thresholds
tic
close all
clear all
clc

SUBGROUP = 0; % If 1, will split the patients by "epilepsy type" - mTLE and not mTLE

load('pat_orig_order.mat')
% pat_ep_type_orig_order: 0 = no mesial SOZs, 1 = pure mesial SOZs, 2 = mesial + other SOZs. 

% Randomize the SOZ labels? ######################
RANDOMIZE = 0; % 1 is TRUE

mA_key = '3mA';

 % Will dictate which data is pulled out of patient structs
% data_type = 'PSD Baseline Z-Scored'; 
% data_type = 'PSD Baseline Mean Subtracted';
data_type = 'RBP Baseline Z-Scored';
% data_type = 'RBP Baseline Mean Subtracted';

% How should patient data be normalized to itself
% intrapat_norm_style = 'No Intra-patient normalization';
% intrapat_norm_style = 'All Z-Score';
intrapat_norm_style = 'Row/Col Z-Score';
% intrapat_norm_style = 'All STD Normalization';
% intrapat_norm_style = 'Row/Col STD Normalized';

source_dir = 'thresholded_dists';
% source_dir = 'X:\000_Data\SPES\results\ccep_psd_distance\gapped_dists\300window_10strideOnlyPT_15gap';
% source_dir = 'X:\000_Data\SPES\results\ccep_psd_distance\gapped_dists\500window_10strideOnlyPT_15gap';

% Get directory contents
contents = dir(source_dir);
mat_idxs = find(contains({contents.name}, '.mat'));
if isempty(mat_idxs); error("%d .MATs found in directory (expected > 0)\n%s",length(mat_idxs), source_dir); end
fprintf("%d MATs found in directory: %s\n",length(mat_idxs), source_dir)
mat_names = string({contents(mat_idxs).name})';
full_mat_paths = fullfile(source_dir,mat_names)';

% Pull out paths for the desired current
full_mat_paths_USED = full_mat_paths(find(contains(mat_names,mA_key)))';


%% Read pat .mats in to prepare for plotting

% Read in one pat to get time windows and freqs used'
S = load(full_mat_paths_USED(1));
windows = S.pat_spes.time_win_ms;
switch data_type
    case {'PSD Baseline Z-Scored','PSD Baseline Mean Subtracted'}
        freqs = S.pat_spes.psd_freqs';
    case {'RBP Baseline Z-Scored', 'RBP Baseline Mean Subtracted'}
        freqs = 1:size(S.pat_spes.rbp_bands,1);
end
dists = S.pat_spes.dist_threshes_mm;
dist_gap = S.dist_gap;

% Initialize plotting variables'
soz_in = nan(length(full_mat_paths_USED),length(dists),size(windows,1),length(freqs));
soz_out = nan(length(full_mat_paths_USED),length(dists),size(windows,1),length(freqs));
pz_in = nan(length(full_mat_paths_USED),length(dists),size(windows,1),length(freqs));
pz_out = nan(length(full_mat_paths_USED),length(dists),size(windows,1),length(freqs));
non_in = nan(length(full_mat_paths_USED),length(dists),size(windows,1),length(freqs));
non_out = nan(length(full_mat_paths_USED),length(dists),size(windows,1),length(freqs));

pats_used = strings(1,length(full_mat_paths_USED));

for i = 1:length(full_mat_paths_USED)
    fprintf("%s\n",full_mat_paths_USED(i))
    
    S = load(full_mat_paths_USED(i));
    curr_struct = S.pat_spes;
    pats_used(i) = curr_struct.patID_clean;
    
    
    % Select the baseline raw mean subtracted OR baseline z-scored data
    % data is ordered (stim x distThresh x time x response x freq)
    switch data_type
        case 'PSD Baseline Z-Scored' 
            data = curr_struct.data_Z;
            data_AllDists = curr_struct.data_Z_AllDists;
            
        case 'PSD Baseline Mean Subtracted'
            data = curr_struct.data_SUB;
            data_AllDists = curr_struct.data_SUB_AllDists;
            
        case 'RBP Baseline Z-Scored'
            data = curr_struct.data_RBP_Z;
            data_AllDists = curr_struct.data_RBP_Z_AllDists;
            
        case 'RBP Baseline Mean Subtracted'
            data = curr_struct.data_RBP_SUB;
            data_AllDists = curr_struct.data_RBP_SUB_AllDists;
    end
    
    stim_soz = curr_struct.stim_soz;
    response_soz = curr_struct.response_soz;
    
    if RANDOMIZE
        stim_soz = stim_soz(randperm(length(stim_soz)));
        response_soz = response_soz(randperm(length(response_soz)));
        disp("WARNING: SOZ LABELS HAVE BEEN RANDOMIZED")
    end
    
    for d = 1:length(dists)
        for w = 1:size(windows,1)
            
            % Use outward normalized data for inward metrics
            % Use inward normalized data for outward metrics
            
            % data is ordered (stim x distThresh x time x response x freq)
            switch intrapat_norm_style

                case 'No Intra-patient normalization'
                    data_inZ = squeeze(data(:,d,w,:,:));
                    data_outZ = squeeze(data(:,d,w,:,:));
                    
                case 'All Z-Score'
                    data_inZ = squeeze((data(:,d,w,:,:) - repmat(mean(data_AllDists(:,1,w,:,:),[1 4],'omitnan'),[size(data_AllDists,1),1,1,size(data_AllDists,4),1])./repmat(std(data_AllDists(:,1,w,:,:),[],[1 4], 'omitnan'),[size(data_AllDists,1),1,1,size(data_AllDists,4),1])));
                    data_outZ = data_inZ;
                    
                case 'Row/Col Z-Score'
                    data_inZ = squeeze((data(:,d,w,:,:) - repmat(mean(data_AllDists(:,1,w,:,:),1,'omitnan'),[size(data_AllDists,1),1,1,1,1]))./repmat(std(data_AllDists(:,1,w,:,:),[],1,'omitnan'),[size(data_AllDists,1),1,1,1,1]));
                    data_outZ = squeeze((data(:,d,w,:,:) - repmat(mean(data_AllDists(:,1,w,:,:),4, 'omitnan'),[1,1,1,size(data_AllDists,4),1]))./repmat(std(data_AllDists(:,1,w,:,:),[],4,'omitnan'),[1,1,1,size(data_AllDists,4),1]));
                    
                case 'All STD Normalization'
                    data_inZ = squeeze(data(:,d,w,:,:)./repmat(std(data_AllDists(:,1,w,:,:),[],[1 4], 'omitnan'),[size(data_AllDists,1),1,1,size(data_AllDists,4),1]));
                    data_outZ = data_inZ;
                    
                case 'Row/Col STD Normalized'
                    data_inZ = squeeze((data(:,d,w,:,:))./repmat(std(data_AllDists(:,1,w,:,:),[],1),[size(data_AllDists,1),1,1,1,1]));
                    data_outZ = squeeze((data(:,d,w,:,:))./repmat(std(data_AllDists(:,1,w,:,:),[],4),[1,1,1,size(data_AllDists,4),1]));
                    
            end
             
            soz_in(i,d,w,:) = squeeze(mean(data_outZ(:,response_soz == 1,:),[1,2],'omitnan'));
            soz_out(i,d,w,:) = squeeze(mean(data_inZ(stim_soz == 1,:,:),[1,2],'omitnan'));
            
            pz_in(i,d,w,:) = squeeze(mean(data_outZ(:,response_soz == 2,:),[1,2],'omitnan'));
            pz_out(i,d,w,:) = squeeze(mean(data_inZ(stim_soz == 2,:,:),[1,2],'omitnan'));
            
            non_in(i,d,w,:) = squeeze(mean(data_outZ(:,response_soz == 3 | response_soz == 0,:),[1,2],'omitnan'));
            non_out(i,d,w,:) = squeeze(mean(data_inZ(stim_soz == 3 | stim_soz == 0,:,:),[1,2],'omitnan'));

        end
    end          
end

toc

%% Agressive power line interference elimination

% Only run this for PSD metrics
soz_in(:,:,:,[57:63,117:123]) = nan;
soz_out(:,:,:,[57:63,117:123]) = nan;
pz_in(:,:,:,[57:63,117:123]) = nan;
pz_out(:,:,:,[57:63,117:123]) = nan;
non_in(:,:,:,[57:63,117:123]) = nan;
non_out(:,:,:,[57:63,117:123]) = nan;

%%  INs & OUTs: Plot metrics per band
close all

edge_types = 'All edge types';
% edge_types = 'No SOZ-SOZ';
% edge_types = 'No SOZ-SOZ, no SOZ-PZ/PZ-SOZ';

band = 1; % band = 8:80; % For PSD this is a range, for RBP this is a single index (e.g. 2 = alpha, refer to 'rbp_bands' in struct)
tws = 1:1;  % for 150 ms windows 1:9 == 5-355 ms
ymin = -0.5;
ymax = 0.5;

distance_for_bar_plot = 1;
jitterAmount = 0.28;

% Distance params
plot_dists = 1:5;
ymin_d = -0.7;
ymax_d = 0.75;

switch edge_types
    
    case 'All edge types'
        % PSD NORMED variables: (patient, dist, tw, freq), values are PSD
        soz_in_dist_band = squeeze(mean(soz_in(:,:,tws,band),[3,4],'omitnan'));
        soz_out_dist_band = squeeze(mean(soz_out(:,:,tws,band),[3,4],'omitnan'));
        pz_in_dist_band = squeeze(mean(pz_in(:,:,tws,band),[3,4],'omitnan'));
        pz_out_dist_band = squeeze(mean(pz_out(:,:,tws,band),[3,4],'omitnan'));
        non_in_dist_band = squeeze(mean(non_in(:,:,tws,band),[3,4],'omitnan'));
        non_out_dist_band = squeeze(mean(non_out(:,:,tws,band),[3,4],'omitnan'));
        
end

soz_color = [0.6350 0.0780 0.1840];
pz_color = [0.8500 0.3250 0.0980];
non_color = [0 0.4470 0.7410];

y_label = sprintf('%s',data_type);

% Bar plots of single distance threshold
if dist_gap > 150
    
    x_label = '';
    
    % INWARD
    f1 = figure(1);
    title_str = sprintf("INWARD %d-%d Hz, %s\n %s, %d-%d ms, Distance > %d mm\n%s\n%s",...
    band(1),band(end), mA_key, data_type, windows(tws(1),1), windows(tws(end),2),dists(distance_for_bar_plot),intrapat_norm_style,edge_types);
    A = soz_in_dist_band;
    B = pz_in_dist_band;
    C = non_in_dist_band;
    f1 = plot_single_dist(f1, A, B, C, distance_for_bar_plot, title_str, x_label, y_label, ymin, ymax, soz_color, pz_color, non_color, jitterAmount);
                           
    % OUTWARD
    f2 = figure(2);
    title_str = sprintf("OUTWARD %d-%d Hz, %s\n %s, %d-%d ms, Distance > %d mm\n%s\n%s",...
    band(1),band(end), mA_key, data_type, windows(tws(1),1), windows(tws(end),2),dists(distance_for_bar_plot),intrapat_norm_style,edge_types);
    A = soz_out_dist_band;
    B = pz_out_dist_band;
    C = non_out_dist_band;
    f2 = plot_single_dist(f2, A, B ,C, distance_for_bar_plot, title_str, x_label, y_label, ymin, ymax, soz_color, pz_color, non_color,jitterAmount);
    
    %  If plotting the patients by epilepsy subtype
    if SUBGROUP
        % Get the patient assignments in proper order
        pat_ep_type_proper_order = get_ep_type_order(mat_names, pat_ep_type_orig_order, pat_orig_order);
        
        mtle_idxs = find(pat_ep_type_proper_order == 1 | pat_ep_type_proper_order == 2);
        nonmtle_idxs = find(pat_ep_type_proper_order == 0);
        
        f3 = figure(3);
        title_str = sprintf("MTLE: INWARD - OUTWARD %d-%d Hz, %s\n %s, %d-%d ms, Distance > %d mm\n%s\n%s",...
        band(1),band(end), mA_key, data_type, windows(tws(1),1), windows(tws(end),2),dists(distance_for_bar_plot),intrapat_norm_style,edge_types);
        A1 = soz_in_dist_band(mtle_idxs) - soz_out_dist_band(mtle_idxs);
        B1 = pz_in_dist_band(mtle_idxs) - pz_out_dist_band(mtle_idxs);
        C1 = non_in_dist_band(mtle_idxs) - non_out_dist_band(mtle_idxs);
        f3 = plot_single_dist(f3, A1, B1, C1, distance_for_bar_plot, title_str, x_label, y_label, -0.4, 0.4, soz_color, pz_color, non_color,jitterAmount);
        
        f4 = figure(4);
        title_str = sprintf("NON-MTLE: INWARD - OUTWARD %d-%d Hz, %s\n %s, %d-%d ms, Distance > %d mm\n%s\n%s",...
        band(1),band(end), mA_key, data_type, windows(tws(1),1), windows(tws(end),2),dists(distance_for_bar_plot),intrapat_norm_style,edge_types);
        A2 = soz_in_dist_band(nonmtle_idxs) - soz_out_dist_band(nonmtle_idxs);
        B2 = pz_in_dist_band(nonmtle_idxs) - pz_out_dist_band(nonmtle_idxs);
        C2 = non_in_dist_band(nonmtle_idxs) - non_out_dist_band(nonmtle_idxs);
        f3 = plot_single_dist(f3, A2, B2, C2, distance_for_bar_plot, title_str, x_label, y_label, -0.4, 0.4, soz_color, pz_color, non_color,jitterAmount);
        

    else
        % INWARD - OUTWARD
        f3 = figure(3);
        title_str = sprintf("INWARD - OUTWARD %d-%d Hz, %s\n %s, %d-%d ms, Distance > %d mm\n%s\n%s",...
        band(1),band(end), mA_key, data_type, windows(tws(1),1), windows(tws(end),2),dists(distance_for_bar_plot),intrapat_norm_style,edge_types);
        A = soz_in_dist_band - soz_out_dist_band;
        B = pz_in_dist_band - pz_out_dist_band;
        C = non_in_dist_band - non_out_dist_band;
        f3 = plot_single_dist(f3, A, B, C, distance_for_bar_plot, title_str, x_label, y_label, -0.4, 0.4, soz_color, pz_color, non_color,jitterAmount);
        
    end
else

    % Plot metrics over distance thresholds

    % INWARD
    f = figure(5);
    clf
    f.Position = [100 100 540 900];
    hold on
    soz_bar_data_IN_allDists = soz_in_dist_band(:,plot_dists);
    pz_bar_data_IN_allDists = pz_in_dist_band(:,plot_dists);
    non_bar_data_IN_allDists = non_in_dist_band(:,plot_dists);
    soz_data_ci95 = 1.96 * nanstd(soz_bar_data_IN_allDists)/sqrt(length(soz_bar_data_IN_allDists));
    pz_data_ci95 = 1.96 * nanstd(pz_bar_data_IN_allDists)/sqrt(length(pz_bar_data_IN_allDists));
    non_data_ci95 = 1.96 * nanstd(non_bar_data_IN_allDists)/sqrt(length(non_bar_data_IN_allDists));
    b=bar(dists(plot_dists),[nanmean(soz_bar_data_IN_allDists); nanmean(pz_bar_data_IN_allDists); nanmean(non_bar_data_IN_allDists)]','FaceColor','flat');
    b(1).CData = soz_color; % SOZ
    b(2).CData = pz_color; % PZ
    b(3).CData = non_color; % Non
    %Dummy objects for legend use
    bh(1) = bar(nan,nan,'FaceColor',soz_color);
    bh(2) = bar(nan,nan,'FaceColor',pz_color);
    bh(3) = bar(nan,nan,'FaceColor',non_color);
    legend(bh,'SOZ','PZ','Non','AutoUpdate','off') 
    % Get the x coordinate of the bars
    ngroups = length(plot_dists);
    nbars = 3;
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    y = [nanmean(soz_bar_data_IN_allDists); nanmean(pz_bar_data_IN_allDists); nanmean(non_bar_data_IN_allDists)];
    err = [soz_data_ci95; pz_data_ci95; non_data_ci95];
    errorbar(x',y',err','LineStyle','none','Color','k','LineWidth',0.5)
    title(sprintf("INWARD %d-%d Hz, %s\n %s, %d-%d ms\n%s\n%s",...
        band(1),band(end), mA_key, data_type, windows(tws(1),1), windows(tws(end),2),intrapat_norm_style,edge_types));
    ylim([ymin_d ymax_d])
    xlabel(sprintf('Distance band (X to X + %d), mm',dist_gap))
    ylabel(y_label)
    hold off
    % % Stats on IN bar chart 
    anova_mat = [soz_bar_data_IN_allDists; pz_bar_data_IN_allDists; non_bar_data_IN_allDists];
    num_pats = size(soz_bar_data_IN_allDists,1);
    anova_mat_filled = fillmissing(anova_mat,'previous',1);
    [p,tbl,stats] = anova2(anova_mat_filled,num_pats,'off');
    annotation('textbox', [0.59, 0.7, 0.1, 0.1], 'String', sprintf("Two-way ANOVA\nDistance:    p=%0.2e\nGroup:        p=%0.2e\nInteraction: p=%0.4f",p(1),p(2),p(3)))

    % OUTWARD
    f = figure(6);
    clf
    f.Position = [100 100 540 900];
    hold on
    soz_bar_data_OUT_allDists = soz_out_dist_band(:,plot_dists);
    pz_bar_data_OUT_allDists = pz_out_dist_band(:,plot_dists);
    non_bar_data_OUT_allDists = non_out_dist_band(:,plot_dists);
    soz_data_ci95 = 1.96 * nanstd(soz_bar_data_OUT_allDists)/sqrt(length(soz_bar_data_OUT_allDists));
    pz_data_ci95 = 1.96 * nanstd(pz_bar_data_OUT_allDists)/sqrt(length(pz_bar_data_OUT_allDists));
    non_data_ci95 = 1.96 * nanstd(non_bar_data_OUT_allDists)/sqrt(length(non_bar_data_OUT_allDists));
    b=bar(dists(plot_dists),[nanmean(soz_bar_data_OUT_allDists); nanmean(pz_bar_data_OUT_allDists); nanmean(non_bar_data_OUT_allDists)]','FaceColor','flat');
    b(1).CData = soz_color; % SOZ
    b(2).CData = pz_color; % PZ
    b(3).CData = non_color; % Non
    %Dummy objects for legend use
    bh(1) = bar(nan,nan,'FaceColor',soz_color);
    bh(2) = bar(nan,nan,'FaceColor',pz_color);
    bh(3) = bar(nan,nan,'FaceColor',non_color);
    legend(bh,'SOZ','PZ','Non','AutoUpdate','off') 
    % Get the x coordinate of the bars
    ngroups = length(plot_dists);
    nbars = 3;
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    y = [nanmean(soz_bar_data_OUT_allDists); nanmean(pz_bar_data_OUT_allDists); nanmean(non_bar_data_OUT_allDists)];
    err = [soz_data_ci95; pz_data_ci95; non_data_ci95];
    errorbar(x',y',err','LineStyle','none','Color','k','LineWidth',0.5)
    title(sprintf("OUTWARD %d-%d Hz, %s\n %s, %d-%d ms\n%s\n%s",...
        band(1),band(end), mA_key, data_type, windows(tws(1),1), windows(tws(end),2),intrapat_norm_style,edge_types));
    ylim([ymin_d ymax_d])
    xlabel(sprintf('Distance band (X to X + %d), mm',dist_gap))
    ylabel(y_label)
    hold off
    % % Stats on IN bar chart 
    anova_mat = [soz_bar_data_OUT_allDists; pz_bar_data_OUT_allDists; non_bar_data_OUT_allDists];
    num_pats = size(soz_bar_data_OUT_allDists,1);
    anova_mat_filled = fillmissing(anova_mat,'previous',1);
    [p,tbl,stats] = anova2(anova_mat_filled,num_pats,'off');
    annotation('textbox', [0.59, 0.7, 0.1, 0.1], 'String', sprintf("Two-way ANOVA\nDistance:    p=%0.2e\nGroup:        p=%0.2e\nInteraction: p=%0.4f",p(1),p(2),p(3)))


    % INWARD - OUTWARD
    f = figure(7);
    clf
    f.Position = [100 100 540 900];
    hold on
    soz_bar_data_INmOUT_allDists = soz_in_dist_band(:,plot_dists) - soz_out_dist_band(:,plot_dists);
    pz_bar_data_INmOUT_allDists = pz_in_dist_band(:,plot_dists) - pz_out_dist_band(:,plot_dists);
    non_bar_data_INmOUT_allDists = non_in_dist_band(:,plot_dists) - non_out_dist_band(:,plot_dists);
    soz_data_ci95 = 1.96 * nanstd(soz_bar_data_INmOUT_allDists)/sqrt(length(soz_bar_data_INmOUT_allDists));
    pz_data_ci95 = 1.96 * nanstd(pz_bar_data_INmOUT_allDists)/sqrt(length(pz_bar_data_INmOUT_allDists));
    non_data_ci95 = 1.96 * nanstd(non_bar_data_INmOUT_allDists)/sqrt(length(non_bar_data_INmOUT_allDists));
    b=bar(dists(plot_dists),[nanmean(soz_bar_data_INmOUT_allDists); nanmean(pz_bar_data_INmOUT_allDists); nanmean(non_bar_data_INmOUT_allDists)]','FaceColor','flat');
    b(1).CData = soz_color; % SOZ
    b(2).CData = pz_color; % PZ
    b(3).CData = non_color; % Non
    %Dummy objects for legend use
    bh(1) = bar(nan,nan,'FaceColor',soz_color);
    bh(2) = bar(nan,nan,'FaceColor',pz_color);
    bh(3) = bar(nan,nan,'FaceColor',non_color);
    legend(bh,'SOZ','PZ','Non','AutoUpdate','off') 
    % Get the x coordinate of the bars
    ngroups = length(plot_dists);
    nbars = 3;
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    y = [nanmean(soz_bar_data_INmOUT_allDists); nanmean(pz_bar_data_INmOUT_allDists); nanmean(non_bar_data_INmOUT_allDists)];
    err = [soz_data_ci95; pz_data_ci95; non_data_ci95];
    errorbar(x',y',err','LineStyle','none','Color','k','LineWidth',0.5)
    title(sprintf("INWARD-OUTWARD %d-%d Hz, %s\n %s, %d-%d ms\n%s\n%s",...
        band(1),band(end), mA_key, data_type, windows(tws(1),1), windows(tws(end),2),intrapat_norm_style,edge_types));
    ylim([-0.5 0.2])
    xlabel(sprintf('Distance band (X to X + %d), mm',dist_gap))
    ylabel(y_label)
    hold off
    % % Stats on IN bar chart 
    anova_mat = [soz_bar_data_INmOUT_allDists; pz_bar_data_INmOUT_allDists; non_bar_data_INmOUT_allDists];
    num_pats = size(soz_bar_data_IN_allDists,1);
    anova_mat_filled = fillmissing(anova_mat,'previous',1);
    [p,tbl,stats] = anova2(anova_mat_filled,num_pats,'off');
    annotation('textbox', [0.59, 0.4, 0.1, 0.1], 'String', sprintf("Two-way ANOVA\nDistance:    p=%0.2e\nGroup:        p=%0.2e\nInteraction: p=%0.4f",p(1),p(2),p(3)))


end