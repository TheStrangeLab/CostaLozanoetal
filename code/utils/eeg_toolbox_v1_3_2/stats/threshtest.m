function [res,t] = threshtest(dat,thresh,duration,lessthan)
%THRESHTEST - Test if data is above/below a thresh for a duration of time.
%
% FUNCTION:
%   [res,t] = threshtest(dat,thresh,duration,lessthan)
%
% INPUT ARGS:
%   dat = p_val;    % Vector of data to test
%   thresh = .001;  % Threshold (test is dat <= thresh)
%   duration = 172; % Duration threshold in samples
%   lessthan = 1;   % (Optional) Direction of test <= or >=
%
% OUTPUT ARGS:
%   res- Binary result of significance
%   t- Vector of start times of significant portions
%

% Check the args
if nargin < 4
  lessthan = 1;
end

% if duration is not a vector, copy to match
duration = duration(:);
if length(duration) ~= size(dat,1)
  duration = ones(size(dat,1),1)*duration(1);
end

% set some defaults
res = zeros(size(dat,1),1);
t = [];

% p-values between 0:1, with higher more significant
if lessthan
  a = (dat <= thresh);
else
  a = (dat >= thresh);
end

% zero pad vector to help find transitions
a = [zeros(size(a,1),1) a zeros(size(a,1),1)];

% find positive and negative inflection points
b = diff(a')';

% locations of the positive inflections
[cposx,cposy] = find(b==1);
cpos = [cposx(:) cposy(:)];
[cnegx,cnegy] = find(b==-1);
cneg = [cnegx(:) cnegy(:)];

if isempty(cpos)
  % no sig points at all
  return
end

% see if must sort rows
if size(dat,1) > 1
  % sort the rows to match start and end
  cpos = sortrows(cpos);
  cneg = sortrows(cneg);
end

% match the durations to the correct rows
durtest = duration(cpos(:,1));

% calc the duration of each epoch
dur = cneg(:,2) - cpos(:,2) - durtest;

% mark the significant ones (>= 0)
res(cpos(dur>=0,1)) = 1;

% save the start and end times of the sig ones
t = [cpos(dur>=0,:) cneg(dur>=0,2)-1];


