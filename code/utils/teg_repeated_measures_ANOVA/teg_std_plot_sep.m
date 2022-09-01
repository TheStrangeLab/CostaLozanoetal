function teg_std_plot_sep(varargin)

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
if length(varargin) > 6,
    levels = varargin{7};
else,
    levels = size(M, 2);
end;
if length(varargin) > 7,
    varnames = varargin{8};
else,
    varnames = {};
end;
if length(varargin) > 8,
    factors = varargin{9};
else,
    factors = {};
end;

if length(factors) >= 3,
    segments = levels(factors(1));
    title00 = varnames{factors(1)};

    leg_levels_id = rec_comb(levels(factors(2:(end - 1))), 0);
    for il = 1:size(leg_levels_id, 1),
        leg0{il} = '';
        for ilf = 1:size(leg_levels_id, 2),
            add0 = [varnames{factors(1 + ilf)} '=' num2str(leg_levels_id(il, ilf))];
            if ilf < (size(leg_levels_id, 2) - 1),
                add0 = [add0 ', '];
            end;
            leg0{il} = [leg0{il} add0];
        end;
    end;
elseif length(factors) == 2,
    segments = 1;
    title00 = '-';

    leg_levels_id = rec_comb(levels(factors(1:(end - 1))), 0);
    for il = 1:size(leg_levels_id, 1),
        leg0{il} = '';
        for ilf = 1:size(leg_levels_id, 2),
            add0 = [varnames{factors(ilf)} '=' num2str(leg_levels_id(il, ilf))];
            if ilf < (size(leg_levels_id, 2) - 1),
                add0 = [add0 ', '];
            end;
            leg0{il} = [leg0{il} add0];
        end;
    end;
else,
    title00 = '';
    segments = 1;
end;

Mbase = M;
ylimext = [Inf -Inf];
for iSeg = 1:segments,
    a = 1 + (iSeg - 1) * size(Mbase, 2) / segments;
    b = a + size(Mbase, 2) / segments - 1;
    M = Mbase(:, a:b);
    
    buf0 = size(M, 2) / 10;
    subplot(1, segments, iSeg);
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
    % t0 = title([title00 title0]);
    t0 = title([title00]);
    set(t0, 'Interpreter', 'none');
    
    yl = ylim;
    ylimext = [min(ylimext(1), yl(1)) max(ylimext(2), yl(2))];
end;

for iSeg = 1:segments,
    subplot(1, segments, iSeg);
    ylim(ylimext);
end;
