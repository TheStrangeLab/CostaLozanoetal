function [F, df1, df2, p, MSM, MSE] = teg_RMF(y_red, nSubj, levels)

y_red = reshape(y_red, nSubj, length(y_red) / nSubj);
df1 = prod(levels - 1);

% Correct for lower-level main effects and interaction.
ADC_prior = [];
ADC_all = [];
for tuple = 1:length(levels),
    tuple_matrix = rec_comb(length(levels) * ones(1, tuple), 1);
    ADC{tuple} = [];
    for iInTuple = 1:size(tuple_matrix, 1),
        ADC0 = get_ADC(y_red, levels, tuple_matrix(iInTuple, :));
        ADC{tuple} = [ADC{tuple} ADC0];
    end;
    if tuple < length(levels),
        ADC_prior = [ADC_prior ADC{tuple}];
    end;
    ADC_all = [ADC_all ADC{tuple}];
end;

% Reduce dimensions via PCA
SADC = cov(ADC_all);
[O2L, EV] = eig(SADC);
ModelX = ADC_all * O2L(:, (end - (df1 - 1)):end);
m = mean(ModelX);
ModelX = ModelX - ones(size(ModelX, 1), 1) * m;

if ~isempty(ADC_prior),
    SADC = cov(ADC_prior);
    [O2L, EV] = eig(SADC);
    f = find(diag(EV) > eps);
    ModelX_prior = ADC_prior * O2L(:, f);
    m = mean(ModelX_prior, 1);
    ModelX_prior = ModelX_prior - ones(size(ModelX_prior, 1), 1) * m;
    % Correct ModelX for ModelX_prior
    y = ModelX;
    X = ModelX_prior;
    b = inv(X'*X)*X'*y;
    pred = X*b;
%     ModelX = ModelX - pred;
else
    ModelX_prior = [];
end;

% Model and test contrasts
[F, df1, df2, p] = teg_test_T2(y_red, ModelX);
MSM = 0;
MSE = 0;

function [lev1, lev2] = inner_parse_diffCode(thisDC)
lev1 = 1;
lev2 = 2;
val = 1;
while val ~= thisDC,
    lev1 = lev1 + 1;
    if lev1 == lev2,
        lev1 = 1;
        lev2 = lev2 + 1;
    end;
    val = val + 1;
end;

function ADC = get_ADC(y_red, levels, factors)

% Create all-differences coder
nDiffCodes = zeros(1, length(levels));
for n = 1:length(factors),
    n0 = levels(factors(n));
    n0 = n0 * (n0 - 1) / 2;
    nDiffCodes(factors(n)) = n0;
end;
nDC2 = nDiffCodes;
f = find(nDC2 == 0);
nDC2(f) = [];
id_matrix = rec_comb(nDC2, 0); % All possible pairwise differences
ADC = [];
for iDiffCode = 1:size(id_matrix, 1),
    diffCode2 = id_matrix(iDiffCode, :);
    diffCode = zeros(1, length(levels));
    diffCode(factors) = diffCode2;
    codeVec = zeros(1, size(y_red, 2));
    prior_factors = 1;
    for iF = 1:length(diffCode),
        if diffCode(iF) == 0,
            prior_factors = prior_factors * levels(iF);
            continue;
        end;
        thisDC = diffCode(iF);
        [lev1, lev2] = inner_parse_diffCode(thisDC);
        tmp = zeros(1, size(y_red, 2));
        tmp = reshape(tmp, prior_factors, length(tmp(:)) / prior_factors);
        prior_factors = prior_factors * levels(iF);
        colsLevel = 1 + mod((1:size(tmp, 2)) - 1, levels(iF));
        tmp = ones(size(tmp, 1), 1) * (colsLevel(:)');
        fPos = find(tmp == lev2);
        fNeg = find(tmp == lev1);
        fZero = find(tmp ~= lev1 & tmp ~= lev2);
        if sum(abs(codeVec)) == 0,
            codeVec(fPos) = 1;
            codeVec(fNeg) = -1;
        else,
            codeVec(fNeg) = -codeVec(fNeg);
            codeVec(fZero) = 0;
        end;
    end;
    ADC = [ADC codeVec(:)];
end;
