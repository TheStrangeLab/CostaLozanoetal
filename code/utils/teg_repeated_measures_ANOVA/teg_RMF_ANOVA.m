function [F, dfM, dfE, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, X_ind] = teg_RMF_ANOVA(expanded0, X0)

nSubj = size(expanded0, 1);
subjM = mean(expanded0, 2);

expanded0_red = zeros(size(expanded0));
for iSubj = 1:size(expanded0, 1),
    y = expanded0(iSubj, :);
    y = y(:) - mean(y(:));
    X_ind = [];
    for iX = 1:size(X0, 2),
        tmp = reshape(X0(:, iX), size(expanded0));
        X_ind = [X_ind tmp(iSubj, :)'];
    end;
    X = X_ind;
    if length(find(X ~= 0)) > 0,
        b = inv(X'*X)*X'*y;
        model = X*b;
    else
        model = zeros(size(y));
    end;
    expanded0_red(iSubj, :) = model;
end;

X = X0;
y = expanded0_red(:);
b = inv(X'*X)*X'*y;
model = X*b;
err = y - model;
SSM = sum(model.^2);
SSE = sum(err.^2);

bModel = b;

try,
    [O2L, ev] = eig(cov(expanded0_red));
catch,
    catcher0 = 1;
end;
fnzev = find(diag(ev) > eps);
if length(fnzev) > 1,
    L = expanded0_red * O2L(:, fnzev);
    S = cov(L);
    eps0 = teg_get_eps(S);
else
    eps0 = 1;
end;

%%
% Girden (1992) recommended that when epsilon is > .75, the Huynh-Feldt 
% correction should be applied and when epsilon is < .75 or nothing is 
% known about sphericity, the Greenhouse-Geisser correction should be 
% applied.
% Girden, E. (1992). ANOVA: Repeated measures. Newbury Park, CA: Sage.
%%

dfM = size(X0, 2);
dfE = (nSubj - 1) * dfM;
dfM_adj = eps0 * dfM;
dfE_adj = eps0 * dfE; 
MSM = SSM / dfM;
MSE = SSE / dfE;
F = MSM / MSE;
p = teg_fsig(F, dfM_adj, dfE_adj);
eta = SSM / (SSM + SSE);
