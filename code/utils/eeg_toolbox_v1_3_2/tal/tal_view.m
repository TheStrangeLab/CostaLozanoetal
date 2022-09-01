function [whichviews,viewpos]=tal_view(drawit,layoutnum)
% TAL_VIEW - Plots the desired views of the brain.
%
% tal_view plots an image of the specified views.
%
% FUNCTION:
%   [whichviews,viewpos]=tal_view(drawit,layoutnum)
%
% INPUT ARGS:
%   drawit - 1 to draw the thing; 0 to just return the positions of the views
%
%   layoutnum - A number, designating the i.d. of the view layout style
%               you want to use. The layouts are hard-coded in this
%               script. They are as follows:
%
%               0 - the horizontal standard view with three surface views
%               1 - the triple-view layout used in the Nature paper
%               2 - a compact layout with all five views
%                   (incl. interhemispherics)
%               3 - six views, like layout 2, but with hipp. view as well
%               4 - four views: rsag, lsag, inf, hipp, in square
%                   arrangement
%

gridcolour=[0.8 0.8 0.8];
forslides=eeganalparams('forslides');

if(nargin<1) drawit=1; end; % draw the thing;
if(nargin<2) layoutnum=0; end; % default layout number--- horizontal,
                               % three views
% First, get all the parameters
[views,viewpos,origins,scales,xcol,ycol,xsign,ysign,XY,hgridlinesx,hgridlinesy,vgridlinesx,vgridlinesy,viewimages]=get_tal_params(layoutnum);
if(isempty(views)) % error-trap this
  fprintf(1,'\nError: Can''t read tal/tal_params.txt\n'); return;
end


switch(layoutnum)
case 0, % the old talmaphorizontal style
%  infsize=2*[0.3282 0.4679];
%  sagsize=2*[0.4332 0.3038];
%  sagshoulder=0.065;
%  viewpos={[0 135],[750 250],[-900 250]};
  whichviews=[1 2 3];
case 1, % the old three-view style from the Nature paper
%  infsize=[0.3282 0.4679];
%  sagsize=[0.4332 0.3038];
%  sagshoulder=0.065;
%  viewpos={[-200 135],[150 1000],[-730 1000]};
  whichviews=[1 2 3];
end 

switch(layoutnum)
 case 0, % the old talmaphorizontal style
  whichviews=[1 2 3];
 case 1, % the old three-view style from the Nature paper
  whichviews=[1 2 3];
 case 2, % horizontal style with the interhemispherics
  whichviews=[1 2 3 4 5];
 case 3, % six views, like layout 2, but with the hippocampal view as well
  whichviews=[1 2 3 4 5 6];
 case 4, % four views: lsag, rsag, inf and hipp, in a square arrangement
  whichviews=[1 2 3 6];
end % switch layoutnum

% get path to pictures
picpath = fileparts(which('tal_view'));

for v=whichviews
  if(drawit)
    if(forslides)
      theview=imread(sprintf('%s/%s_inv.jpg',picpath,viewimages{v}),'jpg');
    else
      theview=imread(sprintf('%s/%s.jpg',picpath,viewimages{v}),'jpg');
    end % reading in the bitmap
    image([viewpos{v}(1) viewpos{v}(1)+size(theview,2)], ...
	  [viewpos{v}(2) viewpos{v}(2)+size(theview,1)],theview);
    hold on;
    % put grid lines and labels
    for i=vgridlinesx{v} % vertical lines
      line([i i],vgridlinesy{v},'Color',gridcolour);
    end
    for i=hgridlinesy{v} % horizontal lines
      line(hgridlinesx{v},[i i],'Color',gridcolour);
    end
  end; % drawit
end % display Tal. axes, cycling through the relevant views

if(drawit)
  AX=axis; thewidth=AX(2)-AX(1); theheight=AX(4)-AX(3); axis off;
  if(forslides) set(gcf,'Color',[0.1 0.1 0.4]); end;
  axis image;axis xy;
  set(gcf,'InvertHardCopy','off','Color',[1 1 1]);
end
