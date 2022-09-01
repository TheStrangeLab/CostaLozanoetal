function [Betw, Betw_labels, Betw_vars_involved] = teg_make_betw_cont(contvar, labels, betw_subj_interaction_depth)

maxTuple = betw_subj_interaction_depth;

[nSubj, nVar] = size(contvar);
for iCol = 1:size(contvar, 2),
    v = contvar(:, iCol);
    u = unique(v);
    if 1 == 0 || length(u) == 2,
        for b = 1:2,
            f = v == u(b);
            contvar(f, iCol) = -1 + 2 * (b - 1);
        end;
    else,
        contvar(:, iCol) = contvar(:, iCol) - mean(contvar(:, iCol));
    end;
end;

Betw = [];
Betw_labels = {};
iBL = 1;
contvar0 = contvar;
Betw_per_tuple = {};
for tuple = 1:min(maxTuple, nVar),
    Betw_per_tuple{tuple} = [];
    id_matrix = rec_comb(nVar * ones(1, tuple), 1);
    for iComb = 1:size(id_matrix, 1),
        newcol = ones(nSubj, 1);
        newlabel = [''];
        bvi = [];
        for iiVar = 1:size(id_matrix, 2),
            iVar = id_matrix(iComb, iiVar);
            bvi = [bvi; iVar];
            newcol = newcol .* contvar(:, iVar);
            newlabel = [newlabel labels{iVar}];
            if iiVar < size(id_matrix, 2),
                newlabel = [newlabel ' x '];
            end;
        end;
        for iprevtup = 1:(tuple - 1),
            % X = [Betw_per_tuple{iprevtup} ones(size(Betw_per_tuple{iprevtup}, 1), 1)];
            X = [Betw_per_tuple{iprevtup}];
            X = X - ones(size(X, 1), 1) * mean(X);
            newcol = newcol - mean(newcol);
            b = inv(X'*X)*X'*newcol;
            pred = X * b;
            newcol = newcol - pred;
        end;
        Betw_per_tuple{tuple} = [Betw_per_tuple{tuple} newcol];
        Betw_labels{iBL} = newlabel;
        Betw_vars_involved{iBL} = bvi;
        iBL = iBL + 1;
    end;
    Betw = [Betw Betw_per_tuple{tuple}];
end;
