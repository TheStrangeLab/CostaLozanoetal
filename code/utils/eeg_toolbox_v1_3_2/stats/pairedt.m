function [sig,t,df] = pairedt(data1, data2)
%PAIREDT Student's paired t-test
%
%

data1 = data1(:);
data2 = data2(:);

if length(data1) ~= length(data2)
  fprintf(1,'PAIREDT: vectors must be the same length.\n');
  return;
end

% paired-samples t-test
cv = cov(data1, data2);
N = size(data1, 1);
df = N-1;

% calculate the t-stat
t = (mean(data1) - mean(data2)) / sqrt( ((var(data1) + var(data2)) - 2*cv(2)) / N );

sig = 1 - tcdf( t, N-1 );
sig = 2 * min(sig, 1-sig);
