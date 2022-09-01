function [ind,peaks] = findpeaks(y)
% FINDPEAKS  Find peaks in real vector.
%  ind = findpeaks(y) finds the indices (ind) which are
%  local maxima in the sequence y.  
%
%  [ind,peaks] = findpeaks(y) returns the value of the peaks at 
%  these locations, i.e. peaks=y(ind);
%
%

% Find the sign of the diffs
%s = sign(diff([0 y(:)' 0]));
s = sign(diff(y(:)'));

% replace the zeros with values to the right
% so that the first val becomes the peak
ind0 = find(s == 0);
ind1 = find(s ~= 0);

% loop over zeros, replacing with first nonzero to the right
for i = ind0
  %r = find(s(i:end)~=0);
  rind = ind1(min(find(ind1 > i)));
  if ~isempty(rind)
    s(i) = s(rind);
  end
end

% get the indices of the maxima
ind=find(diff(s)==-2);

% adjust ind up one
if ~isempty(ind)
  ind = ind + 1;
end

if nargout > 1
  peaks = y(ind);
end


