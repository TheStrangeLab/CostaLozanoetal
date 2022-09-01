function teg_plot_rm_scatter3(contvar_per_sub_per_group, yM, titlestr, xstr, ystr, sepvarlabel)

if size(yM{1}, 2) > 1,
    for n = 1:length(yM),
        yM{n} = mean(yM{n}')';
    end;
    fprintf(['Taking mean of differences in scatter3.\n']);
end;

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[scrsz(3)/10 scrsz(4)/10 0.8*scrsz(3) 0.8*scrsz(4)]);

nPlots = size(yM, 2);
nRows = 1 + floor((nPlots - 1) / 2);
nCols = 1 + mod((nPlots - 1), 2);
if nRows * nCols < size(yM, 2),
    nRows = nRows + 1;
end;
cols = [0 0 0; 1 0 0];
for iCol = 1:length(contvar_per_sub_per_group)
    % subplot(nRows, nCols, iCol);
    cols0 = cols(1 + mod(iCol - 1, 2), :);
    s0 = scatter(contvar_per_sub_per_group{iCol}, yM{iCol}, [], cols0, 'filled');
    hold on;
end;
legend({[sepvarlabel ' low'], [sepvarlabel ' high']});
for iCol = 1:length(contvar_per_sub_per_group)
    % regr line
    cols0 = cols(1 + mod(iCol - 1, 2), :);
    % xlim0 = [min(contvar_per_sub_per_group{iCol}) max(contvar_per_sub_per_group{iCol})];
    xlim0 = xlim;
    [Xplot, regline] = get_regline(contvar_per_sub_per_group{iCol}, yM{iCol}, xlim0);
    l0 = line(Xplot, regline);
    set(l0, 'LineWidth', 2);
    set(l0, 'Color', cols0);
    x0 = xlabel(xstr);
    set(x0, 'Interpreter', 'none');
    y0 = ylabel(ystr);
    set(y0, 'Interpreter', 'none');
end;
t0 = title([titlestr]);
set(t0, 'Interpreter', 'none');
