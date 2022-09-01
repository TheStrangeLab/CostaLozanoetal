function handles=electrodemap(maptype,layoutnum,leads,color,shape,lead_rad)
% ELECTRODEMAP  plot electrode positions on a brain map
%
% handles=electrodemap(maptype,layoutnum,leads,color,shape) returns a
% matrix of handles to the graphics objects painted the rows are each
% view obtained from get_tal_params(layoutnum) 
%
% maptype   : 'num' - lead # labels 
%             'mrkfilled' - plot a filled shape marker for each lead
%             'mrkempty'  - plot an empty shape marker for each lead
%
% layoutnum : the talarack layout to be used for displaying the lead
%             positions 
%
% leads     : a vector of the leads to be displayed
%
% color     : the color with which to draw
%
% shape     : the shape to draw (if maptype is not 'num')
%             known shapes are:
%                'circle'
%                'square'
%                'diamond'
%                'triangle'
%                'itriangle'
%                'star'
%                'cross'
%                'x'
%                'weight'
%                'sideweight'
%
% lead_rad  : the radius of the shape to draw (defaults to 3)
%

% CHANGELOG:
%
% 5/3/2005 - pbs - Fixed plotting bug where more views were getting
% plotted than requested.
%
%



okay=0;

if ~exist('lead_rad','var')
  lead_rad = 3;
end

if(strcmp(maptype,'num')) okay=1; end;
if(strcmp(maptype,'mrkfilled')) okay=1; end;
if(strcmp(maptype,'mrkempty')) okay=1; end;

if(~okay)
   fprintf(1,'\nElectrodemap: I don''t know maptype "%s"\n\n',maptype);
   return;
end

[views,viewpos,origins,scales,xcol,ycol,xsign,ysign,XY]=get_tal_params(layoutnum);

% generate the tal_view
[whichviews,viewpos]=tal_view(0,layoutnum);

% get the coordinates of the electrodes
inferior = [];
lsag = [];
rsag = [];
lsag_int = [];
rsag_int = [];
hipp = [];
for i = 1:length(whichviews)
  switch whichviews(i)
   case 1
    inferior = get_tal_montage('inf');
   case 2
    lsag     = get_tal_montage('lsag');
   case 3
    rsag     = get_tal_montage('rsag');
   case 4
    lsag_int = get_tal_montage('lsag_int');
   case 5
    rsag_int = get_tal_montage('rsag_int');
   case 6
    hipp     = get_tal_montage('hipp');
  end
end

coordinates={inferior,lsag,rsag,lsag_int,rsag_int,hipp};

% Create the electrode values, for subsequent plotting

istemplate=0;

if(strcmp(maptype,'num')) % nothing to do
    func=[leads;ones(1,length(leads))];                                       
elseif(strcmp(maptype,'mrkfilled'))
  func=[leads;ones(1,length(leads))]; % list the leads with a row of ones
                                      % under them
  allcol=ones(length(leads),1)*color; % make the color of the leads
elseif(strcmp(maptype,'mrkempty')) % use white colour for all
  func=[leads;ones(1,length(leads))]; 
  allcol=ones(length(leads),1)*[1 1 1]; % make the color of the leads
  istemplate=1;
end % assembled leads & their colors

% ***** Finally, plot the leads on each of the views *****
handles = ones(length(views),length(leads)).*NaN;

%lead_rad=4; % radius of the lead symbol

for V = 1:length(views) % for each view
  
  if(~isempty(coordinates{V}))
    
    nc=coordinates{V};
    sc=scales(V); 
    O=origins{V}+viewpos{V};
    xi=xcol(V); 
    yi=ycol(V); 
    xfac=xsign(V)*sc; 
    yfac=ysign(V)*sc;
    
    for l=1:size(nc,1) % loop through all lead positions
      
      XXX=O(1)+xfac*nc(l,xi); YYY=O(2)+yfac*nc(l,yi);
      
      if(strcmp(maptype,'num'))
	
	thelead=find(nc(l,1)==func(1,:));
        if(~isempty(thelead))
	  
	  handles(V,thelead)=text(XXX,YYY,sprintf('%i',nc(l,1)));
	  set(handles(V,thelead),'FontWeight','Bold','Color',color,'HorizontalAlignment','center');
      	
	end % if there are leads

      else % not a number map
      
	thelead=find(nc(l,1)==func(1,:));
        if(~isempty(thelead))
	  
          if(istemplate)
            handles(V,thelead) = drawshape(shape,XXX,YYY,lead_rad*sc,[0 0 0]);
          else
            handles(V,thelead) = drawshape(shape,XXX,YYY,lead_rad*sc,allcol(thelead,:));
          end
	  
	end % if there are leads
      end % all the functional maps
    
    end % lead loop
  
  end % if there are leads on this view

end % looping through views


% DRAWSHAPE code

% drawshape(shapename,x,y,r,c,[phi])
% shape names:
% 'circle','square','diamond','triangle','itriangle','star','cross','x','weight','sideweight'
function h = drawshape(shapename,x,y,r,c,phi,edgecolour)

if(strcmp(shapename,'phase'))
  if( (nargin<7) | (isempty(edgecolour)) ) edgecolour=c; end;
else
  if( (nargin<6) | (isempty(phi)) ) edgecolour=c; else edgecolour=phi; end;
end
  
lwidth=0.5;

if(strcmp(shapename,'circle'))  
  t = 0:pi/30:2*pi; X=x+r*sin(t); Y=y+r*cos(t);
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c);  set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'square'))
  X=[x-r x+r x+r x-r x-r];
  Y=[y-r y-r y+r y+r y-r];
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'diamond'))
  a=sqrt(2);
  X=[x-a*r x     x+a*r x     x-a*r];
  Y=[y     y-a*r y     y+a*r y];
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else, h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'triangle'))
  X=[x-r x+r x   x-r];
  Y=[y-r y-r y+r y-r];
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'itriangle'))
  X=[x-r x+r x   x-r];
  Y=[y+r y+r y-r y+r];
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else, h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'star'))
  t=2*pi*[0 2 4 1 3 0]/5; little=1/sqrt(7); R=r*(2/(1+little));
  X=x+R*[0 little*sin(2*pi/10) sin(2*2*pi/10) little*sin(3*2*pi/10) ...
	 sin(4*2*pi/10) little*sin(5*2*pi/10) sin(6*2*pi/10) little* ...
	 sin(7*2*pi/10) sin(8*2*pi/10) little*sin(9*2*pi/10)];
  Y=y+R*[1 little*cos(2*pi/10) cos(2*2*pi/10) little*cos(3*2*pi/10) ...
	 cos(4*2*pi/10) little*cos(5*2*pi/10) cos(6*2*pi/10) little* ...
	 cos(7*2*pi/10) cos(8*2*pi/10) little*cos(9*2*pi/10)];
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'cross'))
  X=x+r*[-5/12 +5/12 +5/12 +1   +1   +5/12 +5/12 -5/12 -5/12 -1   -1   -5/12];
  Y=y+r*[-1   -1   -5/12 -5/12 +5/12 +5/12 +1   +1   +5/12 +5/12 -5/12 -5/12];
  if(c==[0 0 0])
    h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'x'))
  a=r*[-3/12 +3/12 +3/12 +1   +1   +3/12 +3/12 -3/12 -3/12 -1   -1   -3/12];
  b=r*[-1   -1   -3/12 -3/12 +3/12 +3/12 +1   +1   +3/12 +3/12 -3/12 -3/12];
  arot=a*cos(pi/4)+b*sin(pi/4); brot=a*cos(pi/4)-b*sin(pi/4);
  X=x+arot; Y=y+brot;
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'phase'))
  X=[x x+r*cos(phi)];
  Y=[y y+r*sin(phi)];
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',3*lwidth);
  else, set(line(X,Y),'Color',c,'LineWidth',3*lwidth); end
elseif(strcmp(shapename,'weight'))
  a=0.5*r; b=0.3*r;
  X=[x-r x-a x+a x+r x+r x+a x-a x-r x-r];
  Y=[y-r y-b y-b y-r y+r y+b y+b y+r y-r];
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else, h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'sideweight'))
  a=0.5*r; b=0.3*r;
  X=[x-r x+r x+b x+b x+r x-r x-b x-b x-r];
  Y=[y-r y-r y-a y+a y+r y+r y+a y-a y-r];
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else, h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'pentagon'))
  t = 0:pi/2.5:2*pi; X=x+r*sin(t); Y=y+r*cos(t);
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c);  set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'octagon'))
  t = 0:pi/4:2*pi; X=x+r*sin(t); Y=y+r*cos(t);
  if(c==[0 0 0]), h=line(X,Y); set(h,'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c);  set(h,'EdgeColor',edgecolour); end
else fprintf(1,'\nError: %s shapename invalid\n',shapename); end

