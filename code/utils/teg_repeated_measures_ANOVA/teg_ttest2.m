function [t, p, df] = teg_ttest2(X1, X2)

% function [t, p, df] = teg_ttest2(X1, X2)

S = sqrt(var(X1) / length(X1) + var(X2) / length(X2));
t = (mean(X1) - mean(X2)) / S;
df_num = (var(X1) / length(X1) + var(X2) / length(X2)) .^ 2;
df_den = ((var(X1)/length(X1)) .^ 2 / (length(X1) - 1) + (var(X2)/length(X2)) .^ 2 / (length(X2) - 1));
df =  df_num / df_den;

x = (t + sqrt(t.^2 + df)) / (2 * sqrt(t.^2 + df));
z = df / 2;
w = df / 2;
tcdf00 = betainc(x, z, w);
p = 1 - tcdf00;
if p > 0.5,
    p = 1 - p;
end;

