function [p, t, df] = test_diff_corr(contvar1, y_red_M)
% 
% p = 0.666;
% d = 6;
% return;

nIts = 10000;

% y_red_M = y_red_M - mean(y_red_M, 2) * ones(1, 2);
[nSubj, nCol] = size(y_red_M);

cm = corrcoef([contvar1, diff(y_red_M')']);
c = cm(1, 2);

N = nSubj;
t = c / sqrt((1 - c^2) / (N - 2));
df = N - 2;
p = teg_ttest(t, df);

return;

[nSubj, nCol] = size(y_red_M);

cm = corrcoef([contvar1, y_red_M]);
c_w = cm(1, 2:3);
d = c_w(2) - c_w(1);

rand_d = zeros(nIts, 1);
for iIt = 1:nIts,
    rcontvar1 = contvar1(randperm(length(contvar1)));
    cm = corrcoef([rcontvar1, y_red_M]);
    rc_w = cm(1, 2:3);
    rd = rc_w(2) - rc_w(1);
    rand_d(iIt) = rd;
end;
nfapo = length(find(abs(rand_d) > abs(d)));
p = nfapo / nIts;
