function [p, Rbar] = rayleigh_alt(data, alt_n)
% P, Rbar = RAYLEIGH_ALT( DATA, ALT_N )
%
% Returns the p-value for the probability that the data are drawn from a
% uniform distribution rather than a unimodal distribution of unknown mean
% direction. ALT_N refers to an alternate N used to calculate significance.
%
% If DATA is a matrix, returns rayleigh value for each column.

data = data;
n = size(data, 1);
if nargin==2
  n = alt_n;;
end
cols = size(data, 2);

Rbar = []; p = [];
for i=1:cols
  [theta, Rbar(i)] = circmean( data(:, i) );
  Z = n*Rbar(i)^2;
  p(i) = exp(-Z) * (1 + (2*Z - Z^2) / (4*n) - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4) / (288*n^2));
end
