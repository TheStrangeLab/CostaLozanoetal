function [X1, factorStarts1, nColsfactor1, labels, cellsets, factors] = teg_create_ANOVA_dummy(varargin)

levels = varargin{1};
if length(varargin) > 1,
    varnames = varargin{2};
else,
    varnames = cellstr(num2str((1:length(levels))'));
end;

id_matrix = rec_comb(levels, 0);
nComb = size(id_matrix, 1);

X = {};

% Main effects: initial matrix
tuple = 1;
X{tuple} = [];
factorStarts{tuple} = [1];
nColsfactor{tuple} = [];
labels = {};
factors = {};
iLabel = 1;
cellsets = {};
for n = 1:(length(levels) - 1),
    factorStarts{tuple} = [factorStarts{tuple} factorStarts{tuple}(end) + levels(n) - 1];
end;
for iLevelComb = 1:size(id_matrix, 1),
    newrow = [];
    levelComb = id_matrix(iLevelComb, :);
    for iFac = 1:size(levelComb, 2),
        thisLevel = levelComb(iFac) - 1;
        newchunk = zeros(1, levels(iFac) - 1);
        if thisLevel == 0,
            newchunk = -ones(size(newchunk));
        else,
            newchunk(thisLevel) = 1;
        end;
        newrow = [newrow newchunk];
        if iLevelComb == 1,
            nColsfactor{tuple} = [nColsfactor{tuple} length(newchunk)];
        end;
    end;
    X{tuple} = [X{tuple}; newrow];
end;
for iMain = 1:size(id_matrix, 2),
    labels{iLabel} = [varnames{iMain}];
    factors{iLabel} = iMain;
    for iLevel = 1:levels(iMain),
        cellsets{iLabel}{iLevel} = find(id_matrix(:, iMain) == iLevel);
    end;
    iLabel = iLabel + 1;
end;

% Interactions
for tuple = 2:length(levels),
    X{tuple} = [];
    factorStarts{tuple} = [];
    nColsfactor{tuple} = [];
    fs0 = 1;
    tmp0 = length(levels) * ones(1, tuple);
    fac_comb = rec_comb(tmp0, 1);
    for iFacComb = 1:size(fac_comb, 1),
        factorStarts{tuple}(iFacComb) = fs0;
        facs0 = fac_comb(iFacComb, :);
        levels_matrix = rec_comb(levels(facs0) - 1, 0);
        newsubmatrix = [];
        for iLevelComb = 1:size(levels_matrix, 1),
            newvec = ones(nComb, 1);
            for iFactor = 1:size(levels_matrix, 2),
                iCol = factorStarts{1}(facs0(iFactor)) - 1;
                iCol = iCol + levels_matrix(iLevelComb, iFactor);
                newvec = newvec .* X{1}(:, iCol);
            end;
            newsubmatrix = [newsubmatrix newvec];
        end;
        X{tuple} = [X{tuple} newsubmatrix];
        nColsfactor{tuple} = [nColsfactor{tuple} size(newsubmatrix, 2)];
        fs0 = fs0 + size(newsubmatrix, 2);
        labels{iLabel} = []; %[num2str(tuple) '-way '];
        for il = 1:length(facs0),
            labels{iLabel} = [labels{iLabel} '' varnames{facs0(il)}];
            if il < length(facs0),
                labels{iLabel} = [labels{iLabel} ' x '];
            end;
        end;
        factors{iLabel} = facs0;
        
        levels_matrix = rec_comb(levels(facs0), 0);
        for iLevelComb = 1:size(levels_matrix, 1),
            fcs = ones(1, nComb);
            for iFactor = 1:size(levels_matrix, 2),
                f = find(id_matrix(:, facs0(iFactor)) ~= levels_matrix(iLevelComb, iFactor));
                fcs(f) = 0;
            end;
            cellsets{iLabel}{iLevelComb} = find(fcs);
        end;

        iLabel = iLabel + 1;
    end;
end;

X1 = [];
factorStarts1 = [];
nColsfactor1 = [];
fs0 = 0;
for n = 1:length(X),
    X1 = [X1 X{n}];
    factorStarts1 = [factorStarts1 factorStarts{n} + fs0];
    nColsfactor1 = [nColsfactor1 nColsfactor{n}];
    fs0 = fs0 + size(X{n}, 2);
end;
