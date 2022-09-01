% fcircle(shapename,x,y,r,c,[phi])
% shape names:
% 'circle','square','diamond','triangle','itriangle','star','cross','x','weight','sideweight'
function fcircle(shapename,x,y,r,c,phi,edgecolour)


if(strcmp(shapename,'phase'))
  if( (nargin<7) | (isempty(edgecolour)) ) edgecolour=c; end;
else
  if( (nargin<6) | (isempty(phi)) ) edgecolour=c; else edgecolour=phi; end;
end
  
lwidth=0.5;

if(strcmp(shapename,'circle'))
  t = 0:pi/30:2*pi; X=x+r*sin(t); Y=y+r*cos(t);
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c);  set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'square'))
  X=[x-r x+r x+r x-r x-r];
  Y=[y-r y-r y+r y+r y-r];
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'diamond'))
  a=sqrt(2);
  X=[x-a*r x     x+a*r x     x-a*r];
  Y=[y     y-a*r y     y+a*r y];
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else, h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'triangle'))
  X=[x-r x+r x   x-r];
  Y=[y-r y-r y+r y-r];
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'itriangle'))
  X=[x-r x+r x   x-r];
  Y=[y+r y+r y-r y+r];
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else, h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'star'))
  t=2*pi*[0 2 4 1 3 0]/5; little=1/sqrt(7); R=r*(2/(1+little));
  X=x+R*[0 little*sin(2*pi/10) sin(2*2*pi/10) little*sin(3*2*pi/10) ...
	 sin(4*2*pi/10) little*sin(5*2*pi/10) sin(6*2*pi/10) little* ...
	 sin(7*2*pi/10) sin(8*2*pi/10) little*sin(9*2*pi/10)];
  Y=y+R*[1 little*cos(2*pi/10) cos(2*2*pi/10) little*cos(3*2*pi/10) ...
	 cos(4*2*pi/10) little*cos(5*2*pi/10) cos(6*2*pi/10) little* ...
	 cos(7*2*pi/10) cos(8*2*pi/10) little*cos(9*2*pi/10)];
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'cross'))
  X=x+r*[-5/12 +5/12 +5/12 +1   +1   +5/12 +5/12 -5/12 -5/12 -1   -1   -5/12];
  Y=y+r*[-1   -1   -5/12 -5/12 +5/12 +5/12 +1   +1   +5/12 +5/12 -5/12 -5/12];
  if(c==[0 0 0])
    set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'x'))
  a=r*[-3/12 +3/12 +3/12 +1   +1   +3/12 +3/12 -3/12 -3/12 -1   -1   -3/12];
  b=r*[-1   -1   -3/12 -3/12 +3/12 +3/12 +1   +1   +3/12 +3/12 -3/12 -3/12];
  arot=a*cos(pi/4)+b*sin(pi/4); brot=a*cos(pi/4)-b*sin(pi/4);
  X=x+arot; Y=y+brot;
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'phase'))
  X=[x x+r*cos(phi)];
  Y=[y y+r*sin(phi)];
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',3*lwidth);
  else, set(line(X,Y),'Color',c,'LineWidth',3*lwidth); end
elseif(strcmp(shapename,'weight'))
  a=0.5*r; b=0.3*r;
  X=[x-r x-a x+a x+r x+r x+a x-a x-r x-r];
  Y=[y-r y-b y-b y-r y+r y+b y+b y+r y-r];
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else, h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'sideweight'))
  a=0.5*r; b=0.3*r;
  X=[x-r x+r x+b x+b x+r x-r x-b x-b x-r];
  Y=[y-r y-r y-a y+a y+r y+r y+a y-a y-r];
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else, h=patch(X,Y,c); set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'pentagon'))
  t = 0:pi/2.5:2*pi; X=x+r*sin(t); Y=y+r*cos(t);
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c);  set(h,'EdgeColor',edgecolour); end
elseif(strcmp(shapename,'octagon'))
  t = 0:pi/4:2*pi; X=x+r*sin(t); Y=y+r*cos(t);
  if(c==[0 0 0]), set(line(X,Y),'Color',edgecolour,'LineWidth',lwidth);
  else h=patch(X,Y,c);  set(h,'EdgeColor',edgecolour); end
else fprintf(1,'\nError: %s shapename invalid\n',shapename); end
