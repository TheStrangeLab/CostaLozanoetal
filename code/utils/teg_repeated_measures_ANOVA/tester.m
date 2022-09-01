function O = tester

O.falsePos = [];
nIts = 1;
for iIt = 1:nIts,
    levels = [2 3];
    varnames = {};
    varnames2 = {};
    for n = 1:length(levels),
        varnames{n} = ['var' num2str(n)];
    end;
    varnames2{1} = 'var2';
    profile = [-1 0 1 -1 0 1] * 1;
    nSubj = 25;
    nTrials = 200;
    
    M = [];
    M2 = [];
    Betw_label = {'Betw1'};
    Betw = [];
    for iSubj = 1:nSubj,
        depvec = [];
        indepvec = [];
        Betw(iSubj, 1) = randn;
        for iTrial = 1:nTrials,
            indepvec(iTrial, 1) = 1 + floor(rand * levels(1));
            indepvec(iTrial, 2) = 1 + floor(rand * levels(2));
            coder = indepvec(iTrial, 2) + (indepvec(iTrial, 1) - 1) * levels(2);
            depvec(iTrial) = Betw(iSubj) * profile(coder) + randn;
            % depvec(iTrial) = profile(coder) + randn;
        end;
        [M0, dum, N] = rec_combo(depvec(:), indepvec);
        M = [M; [M0 N]];
        %     M = [M; [M0]];
        [M0, dum, N] = rec_combo(depvec(:), indepvec(:, 2));
        M2 = [M2; [M0 N]];
        %     M2 = [M2; [M0]];
    end;
    
    % O.R1 = teg_repeated_measures_ANOVA(M, levels, varnames, Betw, Betw_label, 0, 'test');
    
%     O.R2 = teg_repeated_measures_ANOVA(M2, levels(2), varnames2, Betw, Betw_label, 0, 'test');

    % O.R3 = teg_repeated_measures_ANOVA(M2, levels(2), varnames2, [], {}, 0, 'test');
    O.R3 = teg_repeated_measures_ANOVA(M2, levels(2), varnames2, Betw, Betw_label, 0, 'test');

    O.Betw = Betw;
    O.M = M;
    O.M2 = M2;
    
    f = find(O.R3(:, 4) < 0.05);
    O.falsePos(iIt) = length(f);
end;
