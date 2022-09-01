function [F, df1, df2, p, MSM, MSE, eta] = teg_test_ANOVA(y_red_M, levels)

nSubj = size(y_red_M, 1);
Model = ones(nSubj, 1) * mean(y_red_M);
Error = y_red_M - Model;

SSM = SS(Model(:), mean(Model(:)));
SSE = SS(Error(:), 0);
SST = SSM + SSE;

df1 = prod(levels - 1);
df2 = (nSubj - 1) * df1;
MSM = SSM / df1;
MSE = SSE / df2;
F = MSM / MSE;
p = teg_fsig(F, df1, df2);
eta = SSM / SST;
