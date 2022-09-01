function s = applyToStruct(s,expr,varargin)
%APPLYTOSTRUCT - Perform an operation on each row of a struct.
%
% This function will loop over each row of a structure, using a for
% loop with index ind, which you can use in your expression expr.
% The expression is evaluated after expanding the fieldnames to be
% s(ind).fieldname.
%
% Your expression must contain an assignment to a field of the
% structure for there to be any effect on the returned structure.
%
% FUNCTION:
%   s = applyToStruct(s,expr,varargin);
%
% INPUT ARGS:
%   s = events;  % The structure you want to manipulate
%   expr = 'p = single(p);'; % Expression to expand and evaluate
%   varargin - Variable arguments you can use in the expression
%
% OUTPUT ARGS:
%   s - The struct after evaluating the expression.
%
%

% get the field names
fnames = fieldnames(s);

for f = 1:length(fnames)
  % set the expression to replace
  r_exp = ['\<' fnames{f} '\>'];
  
  % set the replacement
  r_str = ['s(ind).' fnames{f}];
  
  % eval the expression
  expr = regexprep(expr,r_exp,r_str);
end

for ind = 1:length(s)
  eval(expr);
end
