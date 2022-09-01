function y = sem(x,varargin)
% SEM standard error of the mean
% Diego Lozano-Soldevilla DCCN, 21-Mar-2013 12:02:27
if isempty(varargin);
  dim=1;
elseif size(nargin,2)==1;
  dim=varargin{1};
end

n = size(x);
n = n(dim);
y = nanstd(x,[],dim) / sqrt(n);
