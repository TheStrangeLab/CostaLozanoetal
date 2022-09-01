function eps0 = teg_get_eps(S)

if size(S, 1) < 2,
    eps0 = 1;
end;
DS = diag(S);
RS = mean(S, 2);
mS = mean(S(:));

k = size(S, 1);
num = k .^ 2 * (mean(DS(:)) - mS) .^ 2;
den = (k - 1) * (sum(S(:) .^ 2) - 2 * k * sum(RS(:) .^ 2) + k .^ 2 * mS .^ 2);

eps0 = num / den;
