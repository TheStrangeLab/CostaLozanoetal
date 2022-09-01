function signi = teg_fsig(F, df1, df2)

if isnan(F),
    signi = NaN;
    return;
end;

x = df1 * F / (df1 * F + df2);
a = df1 / 2;
b = df2 / 2;
try,
    signi = 1 - betainc(x, a, b);
catch,
    disp('error in teg_fsig');
    signi = 0.666;
end;
signi = 1-fcdf(F, df1, df2);