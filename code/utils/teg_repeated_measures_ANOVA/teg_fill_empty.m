function M = teg_fill_empty(M)

% Empty cells are coded as NaNs.
% Left side are values, right side are counts.
% A single virtual trial is added to replace NaNs.

b = size(M, 2) / 2;
NM = M(:, (b + 1):end);
M = M(:, 1:b);

subjM = zeros(size(M, 1), 1);
subjSD = zeros(size(M, 1), 1);
for ir = 1:size(M, 1),
    row = M(ir, :);
    row(isnan(row)) = [];
    subjM(ir) = mean(row);
    subjSD(ir) = sqrt(var(row(~isnan(row))));
end;
D = M - subjM * ones(1, size(M, 2));
within_eff = zeros(1, size(M, 2));
for ic = 1:size(M, 2),
    col = D(:, ic);
    col(isnan(col)) = [];
    within_eff(ic) = mean(col);
    col = M(:, ic);
    f = find(isnan(col));
    if ~isempty(f),
        M(f, ic) = subjM(f) + within_eff(ic);
        NM(f, ic) = 1;
        fprintf(['Guessing ' num2str(length(f)) ' of ' num2str(size(M, 1)) ' values in column ' num2str(ic), ' row(s) ' num2str(f(:)') '\n']);
    end;
end;

M = [M NM];
