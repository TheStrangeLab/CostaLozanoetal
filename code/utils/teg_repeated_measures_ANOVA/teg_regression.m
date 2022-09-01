function [b, F, df1, df2, p, R2, pred] = teg_regression(varargin)

% function [b, F, df1, df2, p, R2, pred] = teg_regression(X, y, report)

X = varargin{1};

v = var(X);
if isempty(find(v == 0)),
    X = X - ones(size(X, 1), 1) * mean(X);
end;

% [E, V] = eig(X'*X);
% vd = diag(V);
% vd = vd ./ sum(vd);
% f = find(vd > (1e-3 ./ length(vd)));
% X = X * E(:, f);

y = varargin{2};
y = y - mean(y);
if length(varargin) > 2,
    report = varargin{3};
else
    report = '';
end;
b = inv(X' * X) * X' * y;
pred = X * b;
df1 = size(X, 2);
MSM = SS(pred, mean(pred)) / df1;
df2 = (size(y, 1) - 1) - df1;
MSE = SS(pred - y, 0) / df2;
F = MSM / MSE;
p = teg_fsig(F, df1, df2);
R2 = var(pred) / var(y);
if report == 1,
    fprintf(['F(' num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(p) ', R2 = ' num2str(R2)]);
    if p < 0.05,
        fprintf(' ***\n');
    else,
        fprintf(' \n');
    end;
end;
