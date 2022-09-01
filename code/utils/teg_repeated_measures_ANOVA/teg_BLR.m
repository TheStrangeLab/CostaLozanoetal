function Model = teg_BLR(X, y)

% function Model = teg_BLR(X, y)

y = y - min(y);
y = y ./ max(y);
ny = mod(y + 1, 2);
fgr{1} = find(y == 0);
fgr{2} = find(y == 1);

X = X - ones(size(X, 1), 1) * mean(X);
X = X ./ (ones(size(X, 1), 1) * (var(X) .^ (1/2)));

% Fit model
[Model.b_full, Xp] = teg_fit_logistic(X, y);
Model.pred_full = 1 ./ (1 + exp(-Xp * Model.b_full));
acc = length(find(abs(Model.pred_full - y) < abs(Model.pred_full - ny))) / length(y);
Model.acc_full = acc;

% Get pvec
pvec = [];
for ic = 1:size(X, 2),
    [t0, p0, df0] = teg_ttest2(X(fgr{1}, ic), X(fgr{2}, ic));
    pvec(ic) = p0;
end;

% Add variables
[dum, order] = sort(pvec);
incl = [];
pred = mean(y) * ones(size(y));
k = 1;
AICc_crit = getAIC(pred, y, k);
lasterr0 = sum((y - mean(y)) .^ 2);
current_X = [];
for no = 1:length(order),
    fprintf([num2str(no) ' of ' num2str(length(order)) '\n']);
    n = order(no);
    current_X0 = [current_X X(:, n)];
    [b0, X0p] = teg_fit_logistic(current_X0, y);
    pred = 1 ./ (1 + exp(-X0p * b0));
    k = size(X0p, 2);
    AICc = getAIC(pred, y, k);
    err = sum((pred - y) .^ 2);
    if err >= lasterr0,
        continue;
    end;
    [F, pF, df1, df2] = teg_ftest(err, lasterr0, size(current_X0, 2), size(current_X0, 2) - 1, length(y));
    if AICc < AICc_crit & pF < 0.05 / (no + 1),
        lasterr0 = err;
        AICc_crit = AICc;
        current_X = current_X0;
        incl = [incl n];
    end;
end;

X = current_X;

if isempty(X),
    return;
end;

% Re-estimate
[Model.b, Model.Xp] = teg_fit_logistic(X, y);

% Descriptives
Model.variables = incl;
Model.pred = 1 ./ (1 + exp(-Model.Xp * Model.b));
k = size(Model.Xp, 2);
Model.AICc = getAIC(Model.pred, y, k);
acc = length(find(abs(Model.pred - y) < abs(Model.pred - ny))) / length(y);
Model.acc = acc;
Model.R2 = var(Model.pred) / var(y);
[F, p, df1, df2] = teg_ftest(var(Model.pred - y), var(y), size(Model.Xp, 2), 1, length(y));
Model.F = F;
Model.df1 = df1;
Model.df2 = df2;
Model.p = p;

function [b_return, Xp] = teg_fit_logistic(X, y)

% function [b_return, Xp] = teg_fit_logistic(X, y)
% Returns parameter vector b and prediction matrix Xp for use with b.

small_EV = 1e-1;
small_converge = 1e-4;

% Normalize
X = X - ones(size(X, 1), 1) * mean(X);
X = X ./ (ones(size(X, 1), 1) * (var(X) .^ (1/2)));

% Add intercept
X = [ones(size(X, 1), 1) X];

% Save prediction matrix for return.
Xp = X;

% Remove bad trials
vT = var(X');
f = find(vT < eps);
X(f, :) = [];

% Deal with dependent columns if necessary
[O2L, EV] = eig(X' * X);
EV = diag(EV);
onL = 0;
if min(EV) / max(EV) < small_EV,
    fEV = find(EV > max(EV) * small_EV);
    X = X * O2L(:, fEV);
    onL = 1;
end;

% Get initial guess for b from multiple linear regression.
% Solve (X'*X) * b = X' y
A = X' * X;
bb = X' * y;
L = chol(A);
z = inv(L) * bb;
b = inv(L') * z;

err = Inf;
iIt = 0;
maxIt = 1e3;
stopflag = 0;
errvec = [sum((y - mean(y) .^ 2))];
while stopflag == 0,
    iIt = iIt + 1;
    old_b = b;
    old_err = err;
    
    % Make prediction using current b
    pred =  1 ./ (1 + exp(-X * b));
    r = y - pred;
    
    
    % find delta as solution to normal equations:
    % (J'*J) * delta = -J'*J;
    J = [];
    for ir = 1:size(X, 1),
        pi = pred(ir);
        yi = y(ir);
        for ic = 1:size(X, 2),
            xi = X(ir, ic);
            der1 = pi * (1 - pi) * xi;
            J(ir, ic) = der1;
        end;
    end;
    
    % Update b using delta
    % (J'*J) * delta = J'*r;
    A = J' * J;
    bb = J' * r;
    [L, test0] = chol(A);
    if test0 > 0,
        fprintf(['Cholesky factorization failed at iteration ' num2str(iIt) '\n']);
        b = old_b;
        break;
    end;
    z = inv(L) * bb;
    delta = inv(L') * z;
    
    scf = 1;
    b0 = b;
    while 1 == 1,
        b = b0 + scf * delta;
        pred =  1 ./ (1 + exp(-X * b));
        r = y - pred;
        err = sum(r .^ 2);
        if err < old_err,
            break;
        else,
            scf = scf / 2;
        end;
    end;
    errvec = [errvec; err];
    
    % Stopping and housekeeping
    if mean(abs(old_b - b)) / mean(abs(old_b)) < small_converge,
        break;
    end;
    if iIt > 1 & (errvec(end - 1) - err) / errvec(end - 1) < small_converge,
        break;
    end;
    if iIt >= maxIt,
        fprintf(['Max iterations reached.\n']);
        break;
    end;
    old_b = b;
end;
if onL == 1,
    b = O2L(:, fEV) * b;
end;
b_return = b;

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

function AICc = getAIC(fx, y, k)
n = length(y);
sig2 = sum((fx - y).^2) / n;
AIC = 2 * k + n * log(sig2);
AICc = AIC + 2 * k * (k + 1) / (n - k - 1);

function [t, p, df] = teg_ttest2(X1, X2)
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