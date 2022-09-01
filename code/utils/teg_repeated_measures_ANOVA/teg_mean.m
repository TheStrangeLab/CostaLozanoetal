function m = teg_mean(varargin)

M = varargin{1};
if length(varargin) > 1,
    dim0 = varargin{2};
else
    dim0 = 1;
end;

if dim0 == 2,
    M = M';
end;

m = zeros(1, size(M, 2));
for iCol = 1:size(M, 2),
    vec = M(:, iCol);
    f = find(isnan(vec));
    vec(f) = [];
    m(iCol) = mean(vec);
end;

if dim0 == 2,
    m = m';
end;
