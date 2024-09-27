% This is a MATLAB script that performs a series of steps to recreate the 
% illustrations of the manuscript: "Basis profile curve identification to understand 
% electrical stimulation effects in human brain networks"
% by Kai J. Miller, Klaus-Robert Mueller, and Dora Hermes. 
%
% All accompanying files are available at: 
% https://purl.stanford.edu/rc201dv0636.
%
% Please read through the comments at each step of the cell-based sections within 
% the .m-file to understand each step (can execute each cell by using 
% ctrl-enter (PC) or cmd-enter (Mac)). 
% 
% Commentary and notes reference back to the notation of the manuscript.
% 
% Note - This example script assumes your current directory is the bpc_scripts folder.
% kjm 1/2021

%% Add subfolder addon_fxns to path, 
    addpath addon_fxns

%% Load data
    load kjm_ccep_exdata_PHG.mat % this is in the bpc_scripts folder, and is our example data set

%% perform BPC identification 
    [B_struct,excluded_pairs]=bpc_identify(V,pair_types); % this is the function described in the manuscript
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Plot single trials - Make figure for matrix V - as in figure 2D, without stimulation subgroup labels
    load('loc_colormap','cm') % this is a red-blue colormap with gray center (see the "LOC package" tool - kai j miller -- PMID: 17343918 DOI: 10.1016/j.jneumeth.2007.01.019)
    figure, kjm_imagesc_ttb(t_win,1:size(V,2),V)
    % set to gray-centered colormap, and center colorscale on zero
    colormap(cm), set(gca,'clim',max(abs(get(gca,'clim')))*[-1 1]) 
    set(gca,'xtick',[.5:.5:2])
    ylabel('Trial Number'), xlabel('Time from Stimulation (s)')
    colorbar
    set(gcf,'Name','Single stimulation trials (V-matrix)')  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Make figures for CCEPs - matrix G

    %% Calculate group-averaged evoked voltage change matrix, G
    for k=1:length(pair_types)                
        G(k,:)=mean(V(:,pair_types(k).indices),2).';
    end, clear k
    
    %% Plot G matrix
    load loc_colormap
    figure, kjm_imagesc_ttb(t_win,1:size(G,1),G.')
    % set to gray-centered colormap, and center colorscale on zero
    colormap(cm), set(gca,'clim',max(abs(get(gca,'clim')))*[-1 1]) 
    set(gca,'xtick',[.5:.5:2])
    ylabel('Trial Number'), xlabel('Time from Stimulation (s)')
    set(gcf,'Name','Stim-pair averaged matrix (CCEPs // G-matrix)')   
    colorbar
    
    %% Plot individual traces:
    
    % plot example trace - pair 1
        figure, k=1;
        % background lines            
        plot([-.25 2],[0 0],'color',0*[1 1 1],'LineWidth',.5), hold on
        plot([0 0],.9*1000*[-1 1],'color',0*[1 1 1],'LineWidth',.5), hold on  
        % data traces (red on top of white)
        plot(t_win,G(k,:),'w','LineWidth',1.5), hold on
        plot(t_win,G(k,:),'r','LineWidth',1)
        % axes stuff
        set(gca,'ylim',1000*[-1 1],'xlim',[-.25 2])
        box off
        xlabel('Time from Stimulation (s)')
        ylabel('Voltage (\mu V)')
        set(gcf,'Name',['Stim-pair averaged trace plot example (stim-pair ' num2str(k) ')'])  
    
    % plot all example traces
    figure
    for k=1:size(G,1)
        subplot(ceil(size(G,1)/8),8,k)
        % background lines
        plot([-.25 2],[0 0],'color',0*[1 1 1],'LineWidth',.5), hold on
        plot([0 0],.9*1000*[-1 1],'color',0*[1 1 1],'LineWidth',.5), hold on  
        % data traces (red on top of white)
        plot(t_win,G(k,:),'w','LineWidth',1.5), hold on
        plot(t_win,G(k,:),'r','LineWidth',1)
        % axes stuff
        set(gca,'ylim',1000*[-1 1],'xlim',[-.25 2])
        box off
        title(num2str(k))
    end, clear k
    set(gcf,'Name','Stim-pair averaged trace plots (CCEPs)')          
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Make BPC plots

    %% colors for plotting, etc
    maxcols=[[0 1 0];... % green
            [1 1 0];... % yellow
            [1 .7 0];... % orange
            [.7 0 1];... % purple
            [1 0 0];... % red
            [0 0 1];... % blue
            [0 1 1];... % cyan
            ];

    %% plot BPC traces - like in fig 4B, but colored
    traces_fig=figure;
    plot(t_win([1 end]),[0 0],'color',.5*[1 1 1]), hold on
    for q=1:length(B_struct)
        figure(traces_fig)
        hold on, plot(t_win,B_struct(q).curve,'color',maxcols(q,:),'LineWidth',3)
    end, clear q
    set(gca,'YTick',0)
    xlabel('Time from stimulation')
    ylabel('Normalized weight of B.P.C.')
    set(gcf,'Name','Basis Profile Curves')  

    %% calculate projection weights - comment in and out to select different metrics like in supp fig 2
        for bb=1:length(B_struct) % cycle through basis curves
            for n=1:length(B_struct(bb).pairs)   % cycle through pair types represented by this basis curve
                B_struct(bb).plotweights(n)=mean(B_struct(bb).alphas{n}); % mean of alphas
%                 B_struct(bb).plotweights(n)=mean(1-B_struct(bb).ep2{n}./B_struct(bb).V2{n}); % explained variance            
%                 B_struct(bb).plotweights(n)=mean(B_struct(bb).alphas{n}./(B_struct(bb).ep2{n}.^.5)); % alphas normalized by error magnitude
            end, clear n
        end, clear bb    
        
    %% Plot projection weights - like in supp fig S3D, but for PHG site
    figure;
    for q=1:length(B_struct)
        hold on, plot(B_struct(q).pairs,B_struct(q).plotweights,'.','color',maxcols(q,:),'MarkerSize',20)    
    end
    hold on, plot(excluded_pairs,0*excluded_pairs,'k.','MarkerSize',20)
    set(gca,'xlim',[0 length(pair_types)])
    xlabel('Stimulation pair number'), ylabel('Projection weight')    
    set(gcf,'Name','Projection weights')  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Make Brain plots

    gscale=0; %global vs w/i BPC scale (1 is global).
%     gscale=1; %global vs w/i BPC scale (1 is global).

    % viewing angles - comment in and out to select one or the other
    th=280; phi=-20; % viewing angle lateral
%     th=-90; phi=-90; % viewing angle inferior

    % Plot on brain (as in figure 4) 
    B_struct_brainplot(B_struct, brain, th, phi, electrode_locs,interpolated_locs,excluded_pairs,maxcols,gscale)
    set(gcf,'Name','BPC on brain rendering - lateral')  
    

    % viewing angles - comment in and out to select one or the other
%     th=280; phi=-20; % viewing angle lateral
    th=-90; phi=-90; % viewing angle inferior
    
    % Plot on brain (as in figure 4) 
    B_struct_brainplot(B_struct, brain, th, phi, electrode_locs,interpolated_locs,excluded_pairs,maxcols,gscale)
    set(gcf,'Name','BPC on brain rendering - inferior')  
    
    