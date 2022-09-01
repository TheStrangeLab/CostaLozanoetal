function [F, df1, df2, p, MSM, MSE] = teg_RMF(y_red_M, levels)

% y_red_M = reshape(y_red_M, nSubj, length(y_red_M) / nSubj);
nSubj = size(y_red_M, 1);
subjM = mean(y_red_M, 2);
y_red_M = y_red_M - subjM * ones(1, size(y_red_M, 2));

Model = mean(y_red_M);
Model = ones(nSubj, 1) * Model;
Error = y_red_M - Model;

df1 = prod(levels - 1);
df2 = nSubj - df1;

SSM = SS(Model(:));
SSE = SS(Error(:));

MSM = SSM / df1;
MSE = SSE / df2;

F = MSM / MSE;

p = teg_fsig(F, df1, df2);
