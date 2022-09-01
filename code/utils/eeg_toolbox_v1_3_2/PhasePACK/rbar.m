function Rbar = rbar( DATA )
% RBAR = CIRCMEAN( DATA ) 
%
% Function returns the mean resultant length (Rbar) of all of the
% angles in the vector DATA.

if nargin == 0
  help rbar;
  return
end

DATA = DATA(:);
N = length(DATA);

C = sum( cos(DATA) );
S = sum( sin(DATA) );
Rbar = sqrt( C^2 + S^2 )/N;
