function Z = teg_zscore(M)

MeanM = ones(size(M, 1), 1) * mean(M);
SDM = ones(size(M, 1), 1) * sqrt(var(M));
Z = (M - MeanM) ./ SDM;
