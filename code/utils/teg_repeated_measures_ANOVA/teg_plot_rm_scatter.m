function teg_plot_rm_scatter(contvar_per_sub, d, titlestr, xstr, ystr)

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[scrsz(3)/10 scrsz(4)/10 0.6*scrsz(3) 0.6*scrsz(4)]);
s0 = scatter(contvar_per_sub, d, [], [0 0 0], 'filled');
t0 = title(titlestr);
set(t0, 'Interpreter', 'none');
hold on;
% regr line
[Xplot, regline] = get_regline(contvar_per_sub, d, xlim);
l0 = line(Xplot, regline);
set(l0, 'LineWidth', 2);
set(l0, 'Color', [0 0 0]);
x0 = xlabel(xstr);
set(x0, 'Interpreter', 'none');
y0 = ylabel(ystr);
set(y0, 'Interpreter', 'none');
