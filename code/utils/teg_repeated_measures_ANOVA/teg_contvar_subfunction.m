function O = teg_contvar_subfunction(X0, contvar1, cellsets, iPred, M_raw, iBetwCont, O, pCritForFurther, plots, labels, Betw_labels, verbose0, p_fwe, NM, fname)

digprec = 2;


[y_red_M, dum] = teg_inner_recode_raw(M_raw, NM, cellsets, iPred);
M_corr = M_raw;
subjM = mean(M_corr, 2);
M_corr = M_corr - subjM * ones(1, size(M_corr, 2));
[dum, expanded0] = teg_inner_recode_raw(M_corr, NM, cellsets, iPred);
X_red_M = [];
for iX = 1:size(X0, 2),
    [X_red_M0, dum] = teg_inner_recode_raw(reshape(X0(:, iX), size(M_raw)), ones(size(M_raw)), cellsets, iPred);
    X_red_M = [X_red_M X_red_M0(1, :)'];
end;
contvar0 = contvar1{iBetwCont};
contvar0 = contvar0 - mean(contvar0);
contvar0 = contvar0 * ones(1, length(X0) / length(contvar0(:)));
X01 = X0 .* (contvar0(:) * ones(1, size(X0, 2)));
[F, df1, df2, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, X_ind] = teg_RMF_ANOVA(expanded0, X01);

O.R = [O.R; F df1 df2 p MSM MSE];
O.labels{length(O.labels) + 1} = [Betw_labels{iBetwCont} ' x ' labels{iPred}];


if p <= pCritForFurther,
    % Plotting: between, continuous
    fprintf([fname '\t' Betw_labels{iBetwCont} ' x ...\n']);
    prestr = '\t';
    teg_report(prestr, labels, iPred, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, p_fwe, y_red_M, '', verbose0, bModel, X_red_M);
    % Extra for between: correlations with contrast scores
    XX = [];
    for iX = 1:size(X0, 2),
        tmp = X0(:, iX);
        tmp = reshape(tmp, size(expanded0));
        XX = [XX mean(tmp)'];
    end;
%     contrM = expanded0 * XX;
    fprintf([prestr '\tCorrelations:\n']);
    for iCol1 = 1:size(y_red_M, 2),
        for iCol2 = iCol1:size(y_red_M, 2),
            if iCol1 == iCol2,
                vecX = y_red_M(:, iCol1);
                resstr = [prestr '\t\tCell ' num2str(iCol1)];
            else,
                vecX = y_red_M(:, iCol2) - y_red_M(:, iCol1);
                resstr = [prestr '\t\tCell ' num2str(iCol2) ' - ' num2str(iCol1)];
            end;
            vecY = contvar0(:, 1);
            [p, t, df, c] = teg_test_corr([vecX, vecY]);
            resstr = [resstr ', corr = ' num2str(c) ', '];
            resstr = [resstr 't(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p) '\n'];
            if p < 0.05,
                fprintf([resstr]);
            end;
        end;
    end;
    fprintf('\n');
end;
