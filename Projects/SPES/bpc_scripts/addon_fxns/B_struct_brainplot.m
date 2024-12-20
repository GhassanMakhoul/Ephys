function B_struct_brainplot(B_struct, brain, th, phi, electrode_locs,interpolated_locs,excluded_pairs,maxcols,gscale)
% function B_struct_brainplot
% This function is not commented / curated. It goes with the illustration of the methods for 
% "Basis profile curve identification to understand electrical stimulation effects in human brain networks"
% by Kai J. Miller, Klaus-Robert Mueller, and Dora Hermes. 
% Input variables:
%   B_struct: 
%   th, phi: 
%   brain: 
%   electrode_locs:
%   interpolated_locs:
%   excluded_pairs:
%   maxcols:
%   gscale: scale to global maximum across groups (1) or not (any other number)
%   kjm 1/2021

    if exist('gscale')~=1, gscale=0; end % default is to scale each group independently
    
    %%
    figure, ctmr_gauss_plot(brain,[0 0 0],0)
    set(gcf,'color','w')

    %% plot electrodes
    el_size=20;    
    el_add_popout(electrode_locs,'k',el_size,th,phi)
    el_add_popout(electrode_locs,1*[1 1 1],.9*el_size,th,phi)

    %% plot non-significant stim pair sites (in gray)
    el_size=25;    
    el_add_popout(interpolated_locs(excluded_pairs,:),'k',el_size,th,phi)
    el_add_popout(interpolated_locs(excluded_pairs,:),.75*[1 1 1],.9*el_size,th,phi)                  

    %% scaling across groups?
    if gscale==1
        tmp=[];
        for q=1:length(B_struct)
            tmp=[tmp B_struct(q).plotweights];
        end
        scalemax=max(tmp); clear tmp
    end

    %% plot BPCs, colored
    for q=1:length(B_struct)        
        if gscale==1 % scaled to global max 
            color_add_popout_prescaled(interpolated_locs(B_struct(q).pairs,:),B_struct(q).plotweights/scalemax,maxcols(q,:),th,phi)
        else % scaled to BPC max
            color_add_popout_prescaled(interpolated_locs(B_struct(q).pairs,:),B_struct(q).plotweights/max(B_struct(q).plotweights),maxcols(q,:),th,phi)
        end
    end

%% DEPENDENT FUNCTIONS INCLUDED HERE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function color_add_popout_prescaled(locs,wts,maxcol,th,phi)
% function color_add_popout(locs,elcol,msize,th,phi)
% pops electrodes out in th, phi polar coords
% assumes 0<=wts<=1
% all wts should be >=0, if <0, it sets to 0
% kjm 06/20


a_offset=.1*max(abs(locs(:,1)))*[cosd(th-90)*cosd(phi) sind(th-90)*cosd(phi) sind(phi)];

locs=locs+...
    repmat(a_offset,size(locs,1),1);

wts(wts<0)=0;

for k=1:size(locs,1)
        
    k_col=[1 1 1]-([1 1 1]-maxcol)*(wts(k));

            hold on, plot3(locs(k,1),locs(k,2),locs(k,3),'.',...
    'MarkerSize',1.1*(30*abs(wts(k))+20),...
    'Color','k')

    hold on, plot3(locs(k,1),locs(k,2),locs(k,3),'.',...
    'MarkerSize',30*abs(wts(k))+20,...
    'Color',.99*k_col)

end

loc_view(th,phi)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function el_add_popout(locs,elcol,msize,th,phi)
% function el_add_popout(locs,elcol,msize,th,phi)
% pops electrodes out in th, phi polar coords
% kjm 02/11


a_offset=.1*max(abs(locs(:,1)))*[cosd(th-90)*cosd(phi) sind(th-90)*cosd(phi) sind(phi)];

locs=locs+...
    repmat(a_offset,size(locs,1),1);


if ~exist('msize','var')
    msize=8; %marker size
end

if exist('elcol')==0, 
    elcol='r'; %default color if none input
end

hold on, plot3(locs(:,1),locs(:,2),locs(:,3),'.','Color', elcol,'MarkerSize',msize)

loc_view(th,phi)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function loc_view(th, phi)
%function loc_view(theta, phi) 
%this function orients the brain and always puts the lighting behind you
%theta and phi are in degrees, not radians
%make sure the brain plot is your current axis
%this function rotates the brain and lighting to the spherical angle 
%inputted.   it is non-standard b/c of matlab.  so, phi is really "elevation" and not
%phi from standard physics coordinates.  (0,0) is at the back of the brain.  for example: 
%loc_view(180,90) views from the top with the front downward, and
%loc_view(180,90) has the front upward, loc_view(90,0) is from the right,
%maybe i should alter it later so there is the option
%of inputting cartesian, but that can be ambiguous, for example [0 0 1] has ambiguous orientation.
%kjm 08/06


view(th,phi),
%%if you want to go the other way
% [th,phi,r]=cart2sph(view_pt(1),view_pt(1),view_pt(1)); %in radians, but "view" uses degrees with different origin
% th=360*(th+pi/2)/(2*pi);
% phi=360*phi/(pi)

view_pt=[cosd(th-90)*cosd(phi) sind(th-90)*cosd(phi) sind(phi)];
%in order to change the direction of the light:
a=get(gca,'Children');
for i=1:length(a)
    b=get(a(i));
    if strcmp(b.Type,'light') %find the correct child (the one that is the light)
        %object for light is the 2nd element, then use a 
        set(a(i),'Position',view_pt) 
        %or something
    end
end