function [meanvec, idmat, Nvec] = rec_combo(varargin)

depvec = varargin{1};
indepvec = varargin{2};

% Make sure all combinations are included and estimated if necessary
levels = [];
for iVar = 1:size(indepvec, 2),
    u = unique(indepvec(:, iVar));
    levels = [levels length(u)];
end;

indepvec0 = rec_comb(levels, 0);
for iVar = 1:size(indepvec0, 2),
    u = unique(indepvec(:, iVar));
    tmp = indepvec0(:, iVar);
    for iu = 1:length(u),
        f = find(tmp == iu);
        indepvec0(f, iVar) = u(iu);
    end;
end;
depvec0 = NaN * ones(size(indepvec0, 1), 1);
depvec = [depvec0; depvec];
indepvec = [indepvec0; indepvec];
[meanvec, idmat, Nvec] = rec_combo_inner(depvec, indepvec, [], 1, [], [], []);

function [meanvec, idmat, Nvec] = rec_combo_inner(depvec, indepvec, meanvec, depth, idmat, idval, Nvec)

if (depth == size(indepvec, 2)),
    u = unique(indepvec(:, depth));
    for iu = 1:length(u),
        f = find(indepvec(:, depth) == u(iu));
        idval(1, depth) = iu;
        selected = depvec(f);
        selected(find(isnan(selected))) = [];
        z = zscore(selected);
        if length(unique(z)) > 2,
            foutlier = find(abs(z) > 2);
            selected(foutlier) = [];
        end;
        if ~isempty(selected),
            newmeanvec = mean(selected);
            meanvec = [meanvec newmeanvec];
        else,
            meanvec = [meanvec NaN];
        end;
        idmat = [idmat; idval];
        Nvec = [Nvec length(selected)];
    end;
else,
    u = unique(indepvec(:, depth));
    for iu = 1:length(u),
        f = find(indepvec(:, depth) == u(iu));
        idval(1, depth) = iu;
        [meanvec, idmat, Nvec] = rec_combo_inner(depvec(f), indepvec(f, :), meanvec, depth + 1, idmat, idval, Nvec);
    end;
end;
