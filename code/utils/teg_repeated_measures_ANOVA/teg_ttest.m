function [p, t, df] = teg_ttest(varargin)

% function [p, t, df] = teg_ttest(vec)

if length(varargin) == 1,
    vec = varargin{1};
    t = sqrt(length(vec)) * mean(vec) / sqrt(var(vec));
    df = length(vec) - 1;
else
    t = varargin{1};
    df = varargin{2};
end;

p = 0.666;

try,
    x = (t + sqrt(t.^2 + df)) / (2 * sqrt(t.^2 + df));
    z = df / 2;
    w = df / 2;
    tcdf00 = betainc(x, z, w);
    p = 1 - tcdf00;
    if p > 0.5,
        p = 1 - p;
    end;
catch
    return;
end;
