function O = teg_repeated_measures(varargin)

% function O = teg_repeated_measures(M, levels, contvar)
%
% M is observation x variable-combination matrix.
% levels is vector of levels per factor.
%
% Thomas E. Gladwin (2010).

M = varargin{1};
levels = varargin{2};
if length(varargin) > 2,
    contvar = varargin{3};
else
    contvar = [];
end;
[nSubj, nVar] = size(M);

% get id_matrix: contains combination of level per variable (column in M)
id_matrix = rec_comb(levels, 0);

% main contvar effect
if ~isempty(contvar),
    cm = corrcoef([contvar M mean(M, 2)]);
    corrs = cm(2:end, 1);
    [t, df, signi] = corr2t(corrs, size(M, 1));
    for icorr = 1:length(signi),
        if signi(icorr) < 0.05 / length(signi),
            fprintf(['Main contvar effect\n']);
            fprintf(['\t\tContinuous var main effect on level ' num2str(icorr) ', corr: ' num2str(corrs(icorr)) '\n']);
            fprintf(['\t\tt(' num2str(df) ') = ' num2str(t(icorr)) ', p = ' num2str(signi(icorr)) '\n']);
            anysig = 0;
        end;
    end;
end;

% Remove subject effect
subjEffect = mean(M, 2);
M_subj = subjEffect * ones(1, nVar);
M = M - M_subj;

% Find model
explained = zeros(1, size(M, 2));
explained_per_cell = {};
all_profile_vecs = {};
for tuple = 1:length(levels),
    % Find factor-permutations
    tmp0 = length(levels) * ones(1, tuple);
    id_matrix_n = rec_comb(tmp0, 1);
    O.ntuple{tuple} = [];
    
    Current = M - ones(size(M, 1), 1) * explained;
    f_cols_cell = {};
    effectvec = [];
    for iFC = 1:size(id_matrix_n, 1),
        label = num2str(id_matrix_n(iFC, :));
        profile_vecs = [];
        level_comb = rec_comb(levels(id_matrix_n(iFC, :)), 0);
        for iP = 1:size(level_comb, 1),
            selectvec = ones(1, size(M, 2));
            for iiFactor = 1:size(level_comb, 2),
                iFactor = id_matrix_n(iFC, iiFactor);
                f = find(id_matrix(:, iFactor) ~= level_comb(iP, iiFactor));
                selectvec(f) = 0;
            end;
            f_cols_cell{iFC, iP} = find(selectvec);
            mlevel = mean(Current(:, f_cols_cell{iFC, iP}), 2);
            profile_vecs = [profile_vecs mlevel];
        end;
        % Save profile
        all_profile_vecs{tuple, iFC} = profile_vecs;
        profile = mean(profile_vecs);
        O.ntuple{tuple}(iFC, :) = profile;
        effectvec = [effectvec; profile - mean(profile)];
    end;
    % Add effects at this level of combination
    for iFC = 1:size(id_matrix_n, 1),
        explained0 = zeros(1, nVar);
        for iP = 1:size(level_comb, 1),
            f = f_cols_cell{iFC, iP};
            explained0(f) = explained0(f) + effectvec(iFC, iP);
            explained(f) = explained(f) + effectvec(iFC, iP);
        end;
        explained_per_cell{tuple, iFC} = explained0;
    end;
    all_effectvecs{tuple} = effectvec;
end;

for tuple = 1:length(levels),
    % Find factor-permutations
    tmp0 = length(levels) * ones(1, tuple);
    id_matrix_n = rec_comb(tmp0, 1);
    O.ntuple{tuple} = [];
    
    % Test profile_vecs
    for iFC = 1:size(id_matrix_n, 1),
        anysig = 0;
        label = num2str(id_matrix_n(iFC, :));
        profile_vecs = all_profile_vecs{tuple, iFC};
        profile = mean(profile_vecs);
        if length(profile) >= nVar,
            continue;
        end;
        p = prod(levels(id_matrix_n(iFC, :)) - 1);
        k = p + 1;
        
        X = diff(profile_vecs')';
        d = mean(X);
        E = X - ones(nSubj, 1) * d;
        n = nSubj;
        S = cov(X);
        invS = inv(S);
        T2 = n * d(:)' * invS * d(:);
        F = T2 * (n - k + 1) / ((n - 1) * (k - 1));
        df1 = k - 1;
        df2 = n - k + 1;
        signi = teg_fsig(F, df1, df2);
        
        if signi < 0.05,
            fprintf([label '\n']);
            fprintf(['\tProfile: ' num2str(profile) '\n']);
            fprintf(['\tF(' num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(signi) '\n']);
            anysig = 1;
        end;
        
        % contvar
        if ~isempty(contvar),
            cm = corrcoef([contvar profile_vecs]);
            corrs = cm(2:end, 1);
            [t, df, signi] = corr2t(corrs, size(X, 1));
            for icorr = 1:length(signi),
                if signi(icorr) < 0.05 / length(signi),
                    fprintf([label '\n']);
                    fprintf(['\t\tContinuous var main effect on level ' num2str(icorr) ', corr: ' num2str(corrs(icorr)) '\n']);
                    fprintf(['\t\tt(' num2str(df) ') = ' num2str(t(icorr)) ', p = ' num2str(signi(icorr)) '\n']);
                    anysig = 0;
                end;
            end;
        end;
    end;
end;

function id_matrix = rec_comb(levels, perm0)

% id_matrix = rec_comb(levels, perm0)
%
% levels = [n1 n2 ... nm]
% perm0: if 1 only permutations are returned
%
% Thomas Gladwin, 2007.

level_vec = 0 * levels;
id_matrix = rec_loop2_inner(levels, level_vec, [], 1, perm0);
if perm0 == 1,
    frem = [];
    for ir = 1:size(id_matrix, 1),
        u = unique(id_matrix(ir, :));
        if length(u) < size(id_matrix, 2),
            frem = [frem; ir];
        end;
    end;
    id_matrix(frem, :) = [];
end;

function id_matrix = rec_loop2_inner(levels, level_vec, id_matrix, current_factor, perm0)

if current_factor == length(levels),
    loopto = levels(current_factor);
    for iLevel = 1:loopto,
        level_vec(current_factor) = iLevel;
        % check
        beenhere = 0;
        if perm0 == 1,
            s0 = sort(level_vec);
            for iCheck = 1:size(id_matrix, 1),
                s1 = sort(id_matrix(iCheck, :));
                if length(s0) ~= length(s1),
                    continue;
                end;
                if s0 == s1,
                    beenhere = 1;
                    break;
                end;
            end;
        end;
        if beenhere == 0,
            id_matrix = [id_matrix; level_vec];
        end;
    end;
else
    loopto = levels(current_factor);
    for iLevel = 1:loopto,
        level_vec(current_factor) = iLevel;
        id_matrix = rec_loop2_inner(levels, level_vec, id_matrix, current_factor + 1, perm0);
    end;
end;

function [t, df, p] = corr2t(c, N)

t = c ./ sqrt((1 - c.^2) ./ (N - 2));
df = N - 2;
p = 1 - teg_tcdf(t, df);

function signi = teg_fsig(F, df1, df2)

x = df1 * F / (df1 * F + df2);
a = df1 / 2;
b = df2 / 2;
try,
    signi = 1 - betainc(x, a, b);
catch,
    disp('error');
    signi = 0.666;
end;