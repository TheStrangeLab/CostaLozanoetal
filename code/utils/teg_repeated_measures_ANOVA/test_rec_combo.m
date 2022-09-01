function [meanvec, idmat] = test_rec_combo

effects = {};
nEffects = 6;
nLevels = zeros(1, nEffects);
for iEffect = 1:nEffects,
    nLevels(iEffect) = 3; % + floor(2 * rand);
    effects{iEffect} = 1:nLevels(iEffect);
    % effects{iEffect} = (10 ^ (nEffects - iEffect)) * (effects{iEffect} - mean(effects{iEffect}));
    effects{iEffect} = (10 * rand) * (effects{iEffect} - mean(effects{iEffect}));
end;
perf_imperf = 1;
nTrials = 500;
if perf_imperf == 2,
    % Random data
    indepvec = [];
    depvec = [];
    for iTrial = 1:nTrials,
        indeprow = zeros(1, nEffects);
        val = 0;
        for iEffect = 1:nEffects,
            efflev0 = 1 + floor(nLevels(iEffect) * rand);
            eff0 = effects{iEffect}(efflev0);
            indeprow(iEffect) = efflev0;
            val = val + eff0;
        end;
        indepvec = [indepvec; indeprow];
        depvec = [depvec; val];
    end;
else
    % Perfect indepvec matrix: no collinearity
    indepvec = rec_comb(nLevels, 0);
    depvec = [];
    for iTrial = 1:size(indepvec, 1),
        indeprow = indepvec(iTrial, :);
        val = 0;
        for iEffect = 1:nEffects,
            efflev0 = indeprow(iEffect);
            eff0 = effects{iEffect}(efflev0);
            val = val + eff0;
        end;
        depvec = [depvec; val];
    end;
    disp(corrcoef(indepvec));
    % However, if nTrials is set low enough...
    if nTrials < size(indepvec, 1),
        nToRemove = size(indepvec, 1) - nTrials;
        remvec = zeros(size(indepvec, 1), 1);
        remvec(1:nToRemove) = 1;
        remvec = remvec(randperm(length(remvec)));
        indepvec(find(remvec), :) = [];
        depvec(find(remvec), :) = [];
    end;
    disp(corrcoef(indepvec));
end;


testerFac = 1 + floor(rand * nEffects);
% Test all factors
[meanvec, idmat] = rec_combo(depvec, indepvec);
fprintf([num2str(testerFac) ': ' num2str(effects{testerFac}) '\n']);
for iLevel = 1:nLevels(testerFac),
    f = find(idmat(:, testerFac) == iLevel);
    fprintf([num2str(mean(meanvec(f))) ' ']);
end;
fprintf('\n');

% Test only testerFac
[meanvec, idmat] = rec_combo(depvec, indepvec(:, testerFac));
for iLevel = 1:nLevels(testerFac),
    f = find(idmat(:, 1) == iLevel);
    fprintf([num2str(mean(meanvec(f))) ' ']);
end;
fprintf('\n');

% Test testerFac and another factor
otherFac = testerFac + floor((nEffects - 1) * rand);
otherFac = 1 + mod(otherFac, nEffects);
[meanvec, idmat] = rec_combo(depvec, indepvec(:, [testerFac otherFac]));
for iLevel = 1:nLevels(testerFac),
    f = find(idmat(:, 1) == iLevel);
    fprintf([num2str(mean(meanvec(f))) ' ']);
end;
fprintf('\n');
