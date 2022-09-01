function [y_red, pred_red, thisset] = teg_inner_recode_2(y, yN, nSubj, thisset, X0, contvar)

contvar = contvar - mean(contvar);

coder = 0;
for icol = 1:size(X0, 2),
    coder = coder + 10^(icol - 1) * X0(:, icol);
end;
u = unique(coder);
codervec = zeros(size(y));
coderM = [];
for iu = 1:length(u),
    f = find(coder == u(iu));
    codervec(f) = u(iu);
    coderM = [coderM; X0(f(1), :)];
end;

y_red = [];

coder_red = [];
cont_red = [];
for iVarSet = 1:length(thisset),
    y_cell_obs = y(thisset{iVarSet});
    y_cell_obs = reshape(y_cell_obs, nSubj, length(y_cell_obs) / nSubj);
    % Take differing cell counts into account
    y_cell_obsN = yN(thisset{iVarSet});
    y_cell_obsN = reshape(y_cell_obsN, nSubj, length(y_cell_obsN) / nSubj);
    y_cell_obsMN = sum(y_cell_obsN, 2); % Get mean per subject over variable combinations.
    y_cell_obsProp = y_cell_obsN ./ (y_cell_obsMN * ones(1, size(y_cell_obsN, 2)));
    y_cell_obs = y_cell_obs .* y_cell_obsProp;
    tmp = [];
    for iSubj = 1:size(y_cell_obs, 1),
        fNotEmpty = ~isnan(y_cell_obs(iSubj, :));
        tmp(iSubj, 1) = sum(y_cell_obs(iSubj, fNotEmpty), 2);
    end;
    y_cell_obs = tmp;
    y_red = [y_red; y_cell_obs];
    
    coder_red = [coder_red; ones(size(y_cell_obs)) * coderM(iVarSet, :)];
    
    contM = reshape(contvar(thisset{iVarSet}), nSubj, length(thisset{iVarSet}) / nSubj);
    cont_red = [cont_red; contM(:, 1)];
end;

y_red = y_red - mean(y_red);
X = coder_red;
X = X .* (cont_red * ones(1, size(X, 2)));
X = [X ones(size(X, 1), 1)];
b = inv(X'*X)*X'*y_red;
pred_red = X * b;
