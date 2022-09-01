function res = inface(p,v);
%INFACE - Return face point if projection of 3D point p is in face defined by v.
%
%
% FUNCTION:
%   res = inface(p,v);
%
%

% set vertices of face
V1 = v(1,:);
V2 = v(2,:);
V3 = v(3,:);

% make face plane
vec1 = V2-V1;
vec2 = V3-V1;
N = cross(vec1, vec2); % plane is specified by N and one vertex

% see where point p crosses plane on way to origin
t = dot(N, V1)/dot(N, p);
xyz = t*p;

% see if point is within the face
res = [];
if PointInTriangle(xyz,V1,V2,V3)
  res = xyz;
end


% See this site for how to find if point in triangle face
% http://www.blackpawn.com/texts/pointinpoly/default.html

function res  = SameSide(p1,p2, a,b)
cp1 = cross(b-a, p1-a);
cp2 = cross(b-a, p2-a);
res = 0;
if dot(cp1, cp2) >= 0 
  res = 1;
end

function res = PointInTriangle(p, a,b,c)
res = 0;
if SameSide(p,a,b,c) & SameSide(p,b,a,c) & SameSide(p,c,a,b) 
  res = 1;
end


