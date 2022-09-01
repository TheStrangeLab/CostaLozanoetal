function [F, df1, df2, p, MSM, MSE, eta] = teg_RMF_ANOVA(M, NM, factors, levels)

% M = reshape(M, nSubj, length(M) / nSubj);
nSubj = size(M, 1);

% Remove lower-level interactions
% or: correct for all other predictors
oldM = M;
M = inner_remove_lower(M, NM, levels, factors);

% Model and test contrasts
[F, df1, df2, p, MSM, MSE, eta] = teg_test_ANOVA(M, NM, factors, levels);

function M = inner_remove_lower(M, NM, levels, factors)
subjM = mean(M, 2);
M = M - subjM * ones(1, size(M, 2));
columnCodesM = rec_comb(levels, 0);
for tuple = 1:(length(factors) - 1),
% for tuple = 1:length(levels),
    id_matrix = rec_comb(length(levels) * ones(1, tuple), 0);
    for iFacVec = 1:size(id_matrix, 1),
        % get levels of factor combinations
        facVec = id_matrix(iFacVec, :);
        if tuple == length(factors),
            if factors == facVec,
                continue;
            end;
        end;
        combCodes = rec_comb(levels(facVec), 0);
        % Find and remove mean per comb-level
        for iComb = 1:size(combCodes, 1),
            fLevel = 1:size(M, 2);
            for iiFac = 1:size(combCodes, 2),
                iFac = facVec(iiFac);
                f = find(columnCodesM(:, iFac) ~= combCodes(iComb, iiFac));
                fLevel(f) = 0;
            end;
            tmp = zeros(size(M));
            tmp(:, find(fLevel)) = 1;
            ftmp = find(tmp & NM > 0);
            m = M(ftmp) .* NM(ftmp) ./ sum(NM(ftmp));
            M(:, find(fLevel)) = M(:, find(fLevel)) - sum(m);
        end;
    end;
end;
M = M + subjM * ones(1, size(M, 2));

function [F, df1, df2, p, MSM, MSE, eta] = teg_test_ANOVA(M, NM, factors, levels)
Model = zeros(size(M));
IndSc = zeros(size(M));

columnCodesM = rec_comb(levels, 0);
combCodes = rec_comb(levels(factors), 0);
% Find and remove mean per comb-level
for iComb = 1:size(combCodes, 1),
    fLevel = 1:size(M, 2);
    for iiFac = 1:size(combCodes, 2),
        iFac = factors(iiFac);
        f = find(columnCodesM(:, iFac) ~= combCodes(iComb, iiFac));
        fLevel(f) = 0;
    end;
    tmp = zeros(size(M));
    tmp(:, find(fLevel)) = 1;
    ftmp = find(tmp & NM > 0);
    m = M(ftmp) .* NM(ftmp) ./ sum(NM(ftmp));
    Model(:, find(fLevel)) = sum(m);
    IndSc(:, find(fLevel)) = mean(M(:, find(fLevel)), 2) * ones(1, length(find(fLevel)));
end;

Error = IndSc - Model;
subjM = mean(M, 2);
MSubj = subjM * ones(1, size(M, 2));
SSSubj = SS(MSubj(:), mean(M(:)));
SSM = SS(Model(:), mean(M(:)));
SSE_raw = SS(Error(:), 0);
SSE = SSE_raw - SSSubj;
SST = SS(M(:), mean(M(:)));
df1 = prod(levels(factors) - 1);
df2 = (size(M, 1) - 1) * df1;
MSM = SSM / df1;
MSE = SSE / df2;
F = MSM / MSE;
p = teg_fsig(F, df1, df2);
eta = SSM / (SST - SSSubj);
