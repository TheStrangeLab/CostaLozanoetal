function [p, t, df, c] = teg_test_corr(xy)

cm = corrcoef(xy);
c = cm(1, 2);

N = size(xy, 1);
t = c / sqrt((1 - c^2) / (N - 2));
df = N - 2;
p = teg_ttest(t, df);
