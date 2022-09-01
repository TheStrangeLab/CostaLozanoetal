function O = teg_repeated_measures_ANOVA(varargin)

% function O = teg_repeated_measures_ANOVA(M, levels, varnames[, Betw_cont, Betw_labels, plots, fname])
%
% M is observation x variable-combination matrix.
% levels is vector of levels per factor.
% varnames is a cell array of strings.
%
% Other inputs are currently not intended for serious use!
%
% Thomas E. Gladwin (2010).

betw_subj_interaction_depth = 1;
pCritForFurther = 1;
only_show_FDR = 1;
fdr_crit = 0.05;
verbose0 = 1;
plots0 = [];
digprec = 3; % digital precision of printouts.te

M = varargin{1};
levels = varargin{2};
p_fwe = 0.05 / prod(levels);
if size(M, 2) > prod(levels), % Number of observations per cell added
    b = size(M, 2) / 2;
    NM = M(:, (b + 1):end);
    M = M(:, 1:b);
else
    NM = ones(size(M));
end;
[nSubj, nVar] = size(M);
% Remove and store subject effects.
% Take different cell counts into account here.
subjEffect = [];
for iSubj = 1:size(M, 1),
    fNonNaN = ~isnan(M(iSubj, :));
    subjEffect(iSubj, 1) = sum(M(iSubj, fNonNaN) .* NM(iSubj, fNonNaN), 2) ./ sum(NM(iSubj, fNonNaN));
end;
M_subj = subjEffect * ones(1, size(M, 2));%elegant way to do repmat
M_raw = M;
M = M - M_subj;%subtract the mean within participant collapsing all conditions

varnames = varargin{3};
if isempty(varnames),
    varnames = cellstr(num2str((1:length(levels))'));
end;

if length(varargin) > 6,
    fname = varargin{7};
    if exist([fname '.ps']),
        delete([fname '.ps']);
    end;
else
    fname = [];
end;

if length(varargin) > 3,
    contvar = varargin{4};
    contvar0 = contvar;
    Betw_varnames = varargin{5};
    if ~isempty(Betw_varnames),
        [Betw, Betw_labels, Betw_vars_involved] = teg_make_betw_cont(contvar, Betw_varnames, betw_subj_interaction_depth);
        contvar = contvar - ones(size(contvar, 1), 1) * mean(contvar);
        contvar = {};
        % Test between effects of subject scores.
        for iVar = 1:size(Betw, 2),
            Xtmp = [Betw(:, iVar) ones(nSubj, 1)];
            [b, F, df1, df2, p] = teg_regression(Xtmp, subjEffect, 0);
            if p < pCritForFurther,
                fprintf([fname '\t' Betw_labels{iVar} ':']);
                fprintf(['\tF(' num2str(df1, digprec) ', ' num2str(df2, digprec) ') = ' num2str(F, digprec) ', p = ' num2str(p, digprec)]);
                fprintf([' * \n\tb = ' num2str(b(1)) '\t']);
                fprintf('\n');
            end;
            contvar1{iVar} = Betw(:, iVar);
            contvarM = [];
            for iVarM = 1:nVar,
                contvarM = [contvarM Betw(:, iVar)];
            end;
            contvar{iVar} = contvarM(:);
        end;
    else,
        contvar = {};
    end;
else
    contvar = {};
end;

if length(varargin) > 5,
    plots = varargin{6};
else
    plots = 0;
end;
if ~isempty(plots0),
    plots = plots0;
end;

% Get dummy-style matrices.
[X1, factorStarts, nColsFactor, labels, cellsets, factors] = teg_create_ANOVA_dummy(levels, varnames);
y = [];
y_raw = [];
yN = [];
X = [];
for iVar = 1:nVar,
    y = [y; M(:, iVar)];
    y_raw = [y_raw; M_raw(:, iVar)];
    yN = [yN; NM(:, iVar)];
    Xsub = ones(nSubj, 1) * X1(iVar, :);
    X = [X; Xsub];
end;

O.R = [];
O.labels = {};
curr_test = [];

for iPred = 1:length(factorStarts),
    predcols = factorStarts(iPred):(factorStarts(iPred) + nColsFactor(iPred) - 1);
    X0 = X(:, predcols);
    X_red_M = [];
    for iX = 1:size(X0, 2),
        [X_red_M0, dum] = teg_inner_recode_raw(reshape(X0(:, iX), size(M_raw)), ones(size(M_raw)), cellsets, iPred);
        X_red_M = [X_red_M X_red_M0(1, :)'];
    end;
    
    [y_red_M, dum] = teg_inner_recode_raw(M_raw, NM, cellsets, iPred);
    M_corr = M_raw;
    subjM = teg_mean(M_corr, 2);
    M_corr = M_corr - subjM * ones(1, size(M_corr, 2));
    [dum, expanded0] = teg_inner_recode_raw(M_corr, NM, cellsets, iPred);
    [F, df1, df2, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, X_ind] = teg_RMF_ANOVA(expanded0, X0);
    O.R = [O.R; F df1 df2 p MSM MSE eps0];
    O.labels{length(O.labels) + 1} = labels{iPred};
    
    if p <= pCritForFurther,
        teg_report('', labels, iPred, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, p_fwe, y_red_M, fname, verbose0, bModel, X_red_M);
    end;
    
    % Between-continuous
    for iBetwCont = 1:length(contvar),
        if iBetwCont > size(contvar0, 2),
            recursestr = '';
        else,
            recursestr = ' '; % No recurse
        end;
        O = teg_contvar_subfunction(X0, contvar1, cellsets, iPred, ...
            M_raw, iBetwCont, O, pCritForFurther, plots, labels, Betw_labels, verbose0, p_fwe, NM, fname);
        curr_test = curr_test + 1;
    end;
    %close(gcf);
end;

% pvec = O(:, 4);
% [f_fdr, p_fdr] = teg_do_fdr(pvec, fdr_crit);
% fprintf(['FDR p value = ' num2str(p_fdr) ', ' num2str(length(f_fdr)) ' FDR significant test(s) found.\n']);
% 
% if only_show_FDR == 1,
%     pCritForFurther = p_fdr;
% end;