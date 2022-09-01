function x=extend(dat,fs)

% FUNCTION X=EXTEND(DAT,fs)
% returns data set X zeropadded to next multiple of fs
% 
% for matrixes, it will extend each row.
%

if prod(size(dat)) == length(dat)
  % is a vector, so get to correct dim
  dat = dat(:)';
end

% see if we must extend
if(rem(size(dat,2),fs)>0)
  % gonna extend
  newLen = fs*(floor(size(dat,2)/fs)+1);
  x=zeros(size(dat,1),newLen);
  x(:,1:size(dat,2))=dat;
else
  x = dat;
end

