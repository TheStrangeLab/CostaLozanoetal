function teg_std_plot(varargin)

M = varargin{1};
if length(varargin) > 1,
    title0 = varargin{2};
else,
    title0 = 'Figure';
end;
if length(varargin) > 2,
    xlabel0 = varargin{3};
else,
    xlabel0 = 'X';
end;
if length(varargin) > 3,
    ylabel0 = varargin{4};
else,
    ylabel0 = 'Y';
end;
if length(varargin) > 4,
    leg0 = varargin{5};
else,
    leg0 = cellstr(num2str((1:size(M, 2))'));
end;
if length(varargin) > 5,
    ply_sd = varargin{6};
else,
    ply_sd = [];;
end;

buf0 = size(M, 2) / 10;
set(gca, 'NextPlot', 'replacechildren');
colm = [0 0 0; 1 0 0; 0 0 1];
ls = {'-', '--'};
for iLine = 1:size(M, 2),
    l0 = line(1:size(M, 1), M(:, iLine));
    set(l0, 'LineWidth', 2);
    ilm = 1 + mod(iLine - 1, 2);
    ilt = 1 + mod(floor((iLine - 1) / 2), 2);
    set(l0, 'LineStyle', ls{ilt});
    set(l0, 'Color', colm(ilm, :));
end;

xlim([1 - buf0, size(M, 1) + buf0]);
set(gca, 'XTick', 1:size(M, 1));
if size(M, 2) > 1,
    legend(leg0);
    legend(gca, 'Location', 'Best');
end;

for iLine = 1:size(M, 2),
    if ~isempty(ply_sd),
        for iX = 1:size(ply_sd, 1),
            se0 = ply_sd(iX, iLine);
            yl = M(iX, iLine) - se0;
            yh = M(iX, iLine) + se0;
            l0 = line([iX iX], [yl yh]);
            set(l0, 'LineWidth', 2);
            ilm = 1 + mod(iLine - 1, 2);
            ilt = 1 + mod(floor((iLine - 1) / 2), 2);
            set(l0, 'LineStyle', ls{ilt});
            set(l0, 'Color', colm(ilm, :));
        end;
    end;
end;

x0 = xlabel(xlabel0);
set(x0, 'Interpreter', 'none');
y0 = ylabel(ylabel0);
set(y0, 'Interpreter', 'none');
t0 = title(title0);
set(t0, 'Interpreter', 'none');

