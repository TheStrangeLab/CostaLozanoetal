function teg_report(prestr, labels, iPred, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, p_fwe, y_red_M, fname, verbose0, bModel, X_red_M)

% function teg_report(labels, iPred, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, p_fwe, y_red_M, fname, verbose0, bModel, X_red_M)

resstr = [prestr fname ': <' num2str(iPred) '>' '\t' labels{iPred} ': F(' num2str(df1, digprec) ', ' num2str(df2, digprec) ') = ' num2str(F, digprec) ', p = ' num2str(p, digprec) ', eta_p^2 = ' num2str(eta, digprec)];
resstr = [resstr '\n' prestr  '\t\tSSM = ' num2str(SSM, digprec) '\tSSE = ' num2str(SSE, digprec) '\tMSM = ' num2str(MSM, digprec) '\tMSE = ' num2str(MSE, digprec)];
resstr = [resstr '\t epsilon = ' num2str(eps0, digprec)];
fprintf([resstr]);
if p <= p_fwe,
    fprintf(' *** ');
elseif p < 0.05,
    fprintf(' * ');
end;
fprintf('\n');
y_red_Mw = y_red_M - mean(y_red_M, 2) * ones(1, size(y_red_M, 2));
myr = mean(y_red_M);
seyr = var(y_red_Mw) .^ 0.5 / sqrt(size(y_red_M, 1));
fprintf([prestr '\t\ty_cells = [']);
for ib = 1:length(myr),
    fprintf([num2str(myr(ib))]);
    if ib < length(myr),
        fprintf(', ');
    else
        fprintf(']\n');
    end;
end;
fprintf([prestr '\t\tse_cells = [']);
for ib = 1:length(myr),
    fprintf([num2str(seyr(ib))]);
    if ib < length(myr),
        fprintf(', ');
    else
        fprintf(']\n');
    end;
end;
fprintf([prestr '\t\tModel =[\n']);
for ib = 1:size(X_red_M, 2),
    fprintf([prestr '\t\t']);
    fprintf([num2str(X_red_M(:, ib)')]);
    if ib < size(X_red_M, 2),
        fprintf(';\n');
    else
        fprintf(']\n');
    end;
end;
fprintf([prestr '\t\tb = [']);
for ib = 1:length(bModel),
    fprintf([num2str(bModel(ib))]);
    if ib < length(bModel),
        fprintf(', ');
    else
        fprintf(']\n');
    end;
end;
if verbose0 > 0,
    if length(myr) > 2,
        % Paired t-tests
        fprintf([prestr '\tPost-hoc differences:\n']);
        for col1 = 1:size(y_red_M, 2),
            for col2 = 1:(col1 - 1),
                vec = y_red_M(:, col1) - y_red_M(:, col2);
                [p, t, df] = teg_ttest(vec);
                if p < 0.05,
                    tresstr = [num2str(col1) ' - ' num2str(col2) ': t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p)];
                    fprintf([prestr '\t\t' tresstr '\n']);
                end;
            end;
        end;
        % vs 0
        fprintf('\tPost-hoc test vs 0:\n');
        for col1 = 1:size(y_red_M, 2),
            vec = y_red_Mw(:, col1);
            [p, t, df] = teg_ttest(vec);
            if p < 0.05,
                tresstr = [num2str(col1) ': t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p)];
                fprintf([prestr '\t\t' tresstr '\n']);
            end;
        end;
    else,
        fprintf('\n');
    end;
    fprintf('\n');
else,
    fprintf('\n');
end;
