function [F, df1, df2, p, MSM, MSE, eta] = teg_test_T2(y_red, ModelX)

ContrastsM = y_red * ModelX;
k = size(ContrastsM, 2) + 1;
meandiffvec = mean(ContrastsM)';
S = cov(ContrastsM);
n = size(y_red, 1);
T2 = n * meandiffvec' * inv(S) * meandiffvec;
F = T2 * (n - k + 1) / ((n - 1) * (k - 1));
df1 = k - 1;
df2 = n - df1;
p = teg_fsig(F, df1, df2);

MSM = T2;
MSE = det(S);
eta = df1 * F / (df1 * F + df2);

doorstop = 1;
