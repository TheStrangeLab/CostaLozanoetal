function [F, p, df1, df2] = teg_ftest(Err_new, Err_old, p_new, p_old, n)

% function [F, p, df1, df2] = teg_ftest(Err_new, Err_old, p_new, p_old, n)
%
% F test for significance of the error variance (Err) decrease from 
% a p_old parameter model to a p_new parameter model.
% Include intercept when counting parameters.
% n is number of observations.

RSS_old = Err_old * (n - 1);
RSS_new = Err_new * (n - 1);

df1 = p_new - p_old;
df2 = n - p_new;

F = ((RSS_old - RSS_new) / df1) / (RSS_new / df2);

x = df1 * F / (df1 * F + df2);
a = df1 / 2;
b = df2 / 2;
p = 1 - betainc(x, a, b);
