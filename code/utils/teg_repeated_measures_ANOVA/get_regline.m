function [Xplot, regline] = get_regline(contvar_per_sub, d, xlim0)
X = [contvar_per_sub ones(size(contvar_per_sub, 1), 1)];
y = d;
b = inv(X'*X)*X'*y;
Xplot = [xlim0(:) [1; 1]];
regline = Xplot * b;
Xplot = Xplot(:, 1);
