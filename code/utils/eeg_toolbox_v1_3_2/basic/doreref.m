function doreref(r)
%DOREREF - Wrapper function for reref
%
% Not integrated into toolbox.
%
% FUNCTION doreref(r)
%
% Wrapper to the reref script, making it easier to define the grids
% of electrodes.
%
% R is made up of multiple rows of two integers each, indicating the
% beginning and end of each grid.
%
% INPUT ARGS:
%   r = [1,8;
%       9,24;
%       25,40;
%       41,48;
%       49,56;
%       57,64;];
%
% OUTPUT:
%   Rereferenced electrode data will be put in the dat/dat
%   directory.
%
% REQUIREMENTS:
%   - 'tal/good_leads.txt'
%   - 'tal/leads.txt'
%


% set the vars and read in good leads
al = getleads('tal/good_leads.txt');
lfile = 'tal/leads.txt';
outdir = 'dat/';

% get the weights
w = ones(size(al));
for i = 1:size(r,1)
  % find the index of those good leads which are in this grid
  idx = find(al>=r(i,1) & al<=r(i,2));
  % weight them appropriately
  w(idx) = w(idx)/length(idx);
end

% combine the good leads and weights
al = [al;w];

% reref, yo
reref(outdir, lfile, al);


