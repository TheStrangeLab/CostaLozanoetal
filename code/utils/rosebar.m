function h = rosebar(theta,x,varargin)
edges = sort(rem([(x(2:end)+x(1:end-1))/2 (x(end)+x(1)+2*pi)/2],2*pi));
edges = [edges edges(1)+2*pi];

nn = theta(:);
[m,n] = size(nn);
mm = 4*m;
r = zeros(mm,n);
r(2:4:mm,:) = nn;
r(3:4:mm,:) = nn;

t = zeros(mm,1);
t(2:4:mm) = edges(1:m);
t(3:4:mm) = edges(2:m+1);

h = polar(t,r);
