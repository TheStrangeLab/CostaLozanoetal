% calculate standard error
% 
% [s] = sterror(x)
% 
% output s: standard error of x
% input  x: vector of scalar values

 
function [s] = sterror(x);

 
s = sqrt(var(x))./sqrt(length(x));

% Stephan