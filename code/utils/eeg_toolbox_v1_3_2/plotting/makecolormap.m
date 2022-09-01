function map = makecolormap(startcol,endcol,numpoints)
%MAKECOLORMAP - Generate a gradient between two RGB colors.
%
% FUNCTION:
%   map = makecolormap(startcol,endcol,numpoints)
%
%
%
%





% get the delta values to use
delta = (endcol-startcol)./(numpoints-.99999);

% handle no change
delta(delta==0) = 1;

map = [(startcol(1):delta(1):endcol(1)).*ones(1,numpoints);
       (startcol(2):delta(2):endcol(2)).*ones(1,numpoints);
       (startcol(3):delta(3):endcol(3)).*ones(1,numpoints);]';

% matlab rounding weirdness
%x = 1000;
%map = round(map.*x)./x;
