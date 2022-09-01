M = [2 4 6 10 8 6];
M = ones(40, 1) * M;
M = M + randn(size(M))
levels = [2 3];
varnames = {'Fac1', 'Fac2'};
O = teg_repeated_measures_ANOVA(M, levels, varnames)



levels = [3 2 2];
varnames = {'drug', 'hemi', 'ipsicont'};
O = teg_repeated_measures_ANOVA(SPSS, levels, varnames)