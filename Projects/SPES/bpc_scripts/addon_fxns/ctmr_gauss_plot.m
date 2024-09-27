function ctmr_gauss_plot(cortex,electrodes,weights)
% function [electrodes]=ctmr_gauss_plot(cortex,electrodes,weights)
% projects electrode locations onto their cortical spots in the 
% left hemisphere and plots about them using a gaussian kernel
% for only cortex use: 
% ctmr_gauss_plot(cortex,[0 0 0],0)
% rel_dir=which('loc_plot');
% rel_dir((length(rel_dir)-10):length(rel_dir))=[];
% addpath(rel_dir)
%   Created by:
%   K.J. Miller & D. Hermes, 
%   Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%   Version 1.1.0, released 26-11-2009

%%
    %load in colormap
    load('loc_colormap')

    brain=cortex.vert;
    v='l';

    if length(weights)~=length(electrodes(:,1))
        error('you sent a different number of weights than electrodes (perhaps a whole matrix instead of vector)')
    end
%% gaussian "cortical" spreading parameter - in mm, so if set at 10, its 1 cm
%- distance between adjacent electrodes
    gsp=50;

    c=zeros(length(cortex(:,1)),1);
    for i=1:length(electrodes(:,1))
        b_z=abs(brain(:,3)-electrodes(i,3));
        b_y=abs(brain(:,2)-electrodes(i,2));
        b_x=abs(brain(:,1)-electrodes(i,1));
    %     d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2).^.5)/gsp^.5); %exponential fall off 
        d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian 
        c=c+d';
    end

%%
% c=(c/max(c));
a=tripatch(cortex, 'nofigure', c');
shading interp;
a=get(gca);
%%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
%%RELEVANCE OF CHANNEL FOR COMPARISON's ACROSS CORTICES
d=a.CLim;
set(gca,'CLim',[-max(abs(d)) max(abs(d))])
l=light;
colormap(cm)
lighting gouraud; %play with lighting...
% material dull;
material([.3 .8 .1 10 1]);
axis off
set(gcf,'Renderer', 'zbuffer')

if v=='l'
view(270, 0);
set(l,'Position',[-1 0 1])        
elseif v=='r'
view(90, 0);
set(l,'Position',[1 0 1])        
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function handle=tripatch(struct, nofigure, varargin)
% TRIPATCH handle=tripatch(struct, nofigure)
%   Created by: ???

if nargin<2 | isempty(nofigure)
   figure
end% if
if nargin<3
   handle=trisurf(struct.tri, struct.vert(:, 1), struct.vert(:, 2), struct.vert(:, 3));
else
   if isnumeric(varargin{1})
      col=varargin{1};
      varargin(1)=[];
      if [1 3]==sort(size(col))
         col=repmat(col(:)', [size(struct.vert, 1) 1]);
      end% if
      handle=trisurf(struct.tri, struct.vert(:, 1), struct.vert(:, 2), struct.vert(:, 3), ...
         'FaceVertexCData', col, varargin{:});
      if length(col)==size(struct.vert, 1)
         set(handle, 'FaceColor', 'interp');
      end% if
   else
      handle=trisurf(struct.tri, struct.vert(:, 1), struct.vert(:, 2), struct.vert(:, 3), varargin{:});
   end% if
end% if
axis tight
axis equal
hold on
if version('-release')>=12
   cameratoolbar('setmode', 'orbit')
else
   rotate3d on
end% if    
