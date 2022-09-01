function teg_plot_rm_scatter2(contvar_per_sub, yM, titlestr, xstr, ystr)

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[scrsz(3)/10 scrsz(4)/10 0.8*scrsz(3) 0.8*scrsz(4)]);

nPlots = size(yM, 2);
nRows = 1 + floor((nPlots - 1) / 2);
nCols = 1 + mod((nPlots - 1), 2);
if nRows * nCols < size(yM, 2),
    nRows = nRows + 1;
end;
cols = [0 0 0; 0 0 0];
for iCol = 1:size(yM, 2)
    subplot(nRows, nCols, iCol);
    cols0 = cols(1 + mod(iCol - 1, 2), :);
    s0 = scatter(contvar_per_sub, yM(:, iCol), [], cols0, 'filled');
    t0 = title([titlestr ' ' num2str(iCol)]);
    set(t0, 'Interpreter', 'none');
    hold on;
    % regr line
    [Xplot, regline] = get_regline(contvar_per_sub, yM(:, iCol), xlim);
    l0 = line(Xplot, regline);
    set(l0, 'LineWidth', 2);
    set(l0, 'Color', cols0);
    x0 = xlabel(xstr);
    set(x0, 'Interpreter', 'none');
    y0 = ylabel(ystr);
    set(y0, 'Interpreter', 'none');
end;
