function f = plot_rosegroup_pacoi(dt)

myFontSize  = 0.4;

% f = figure;%('units', 'centimeters', 'position', [5, 5, 8.7, 12], 'visible', 'on')
% axis off;

f = figure('units', 'centimeters', 'position', dt.figEx);
% if strcmp(dt.visible, 'off')
%     set(f, 'visible', 'off');
% end
axis off;

%----- angle-circle
ax1 = axes('units', 'centimeters', 'position', dt.boxEx);

hold on;

% plot angle lines
for iAng = 0:30:165
    plot([-cosd(iAng), cosd(iAng)], [-sind(iAng), sind(iAng)], ':', ...
        'Color', rgb('gray'), 'LineWidth', 1);
end

% plot circle
tmpx = -pi:0.001:pi;
tmpc = linspace(0, 1, length(tmpx));
patch([cos(tmpx), 0.95 .* fliplr(cos(tmpx))], ...
    [sin(tmpx), 0.95 .* fliplr(sin(tmpx))], [tmpc, fliplr(tmpc)], ...
    'EdgeColor', 'interp');
% plot angle text
t1 = text(1.05, 0, dt.angLabels{1}, ...
    'HorizontalAlignment', 'left');
t2 = text(0, 1.125, dt.angLabels{2}, ...
    'HorizontalAlignment', 'center');
t3 = text(-1.05, 0, dt.angLabels{3}, ...
    'HorizontalAlignment', 'right');
t4 = text(0, -1.125, dt.angLabels{4}, ...
    'HorizontalAlignment', 'center');
% indicate max pac value
t5 = text(0.175, 0.15, ['max ', num2str(max(dt.FR), '%0.2f'), '/', num2str(max(dt.FR1), '%0.2f')], ...
     'units', 'normalized', ...
     'horizontalalignment', 'right', 'verticalalignment', 'top');
% % indicate number of subjects used for the analysis
 t6 = text(1.05, 0.15, 'n=9', ...
         'units', 'normalized', ...
         'horizontalalignment', 'left', 'verticalalignment', 'top');
set([t1, t2, t3, t4, t5, t6], ...0`l
    'FontUnits', 'centimeters', 'FontSize', myFontSize);
hold off;
set(gca, ...
    'xlim', [-1.05, 1.05], 'ylim', [-1.05, 1.05]);
% flip y-axis because the arena is programmed such that negative yaw-values
% point north
if dt.bReverseY == true
    set(gca, ...
        'ydir', 'reverse');
end
myCM = circshift(hsv, size(hsv, 1) / 2); % so that east (0°) is red
colormap(ax1, myCM);
axis off;
%----- maximum pacoi value
ax2 = axes('units', 'centimeters', 'position', ax1.Position);
maxFR = max(dt.FR);

hold on;
for iA = 1:size(dt.angleCenters, 1)
    thisAngles  = transpose([0, dt.angleCenters(iA) - deg2rad(dt.angularRes / 2):0.001:dt.angleCenters(iA) + deg2rad(dt.angularRes / 2), 0]);
    thisRadii   = [0; repmat(dt.FR(iA), numel(thisAngles) - 2, 1); 0] ./ maxFR; % normalize in relation to the maximum pacoi value
    [x, y]      = pol2cart(thisAngles, thisRadii);
    patch(x, y, [1, 0, 0], ...
        'FaceAlpha', '0.3', 'EdgeColor', 'none');
    if dt.angularRes < 15
        set(p1, 'LineWidth', 0.5);
    end
end
% circular mean
circM               = circ_mean(dt.angleCenters, dt.FR);
[circMx, circMy]    = pol2cart(circM, 1);
plot([0, circMx], [0, circMy], '-', ...
    'linewidth', 3, 'color', [1, 0, 0]);
hold off;
set(gca, ...
    'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
if dt.bReverseY == true
    set(gca, ...
        'ydir', 'reverse'); % flip y-axis
end
axis off;
hold on
maxFR1 = max(dt.FR1);
hold on
for iA = 1:size(dt.angleCenters, 1)
    thisAngles  = transpose([0, dt.angleCenters(iA) - deg2rad(dt.angularRes / 2):0.001:dt.angleCenters(iA) + deg2rad(dt.angularRes / 2), 0]);
    thisRadii   = [0; repmat(dt.FR1(iA), numel(thisAngles) - 2, 1); 0] ./ maxFR1; % normalize in relation to the maximum pacoi value
    [x, y]      = pol2cart(thisAngles, thisRadii);

    patch(x, y, [0, 0, 0], ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    if dt.angularRes < 15
        set(p1, 'LineWidth', 0.5);
    end
end
circM1               = circ_mean(dt.angleCenters, dt.FR1);
[circMx1, circMy1]    = pol2cart(circM1, 1);
plot([0, circMx1], [0, circMy1], '-', ...
    'linewidth', 3, 'color', [0, 0, 0]);
hold off;
set(gca, ...
    'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2]);
if dt.bReverseY == true
    set(gca, ...
        'ydir', 'reverse'); % flip y-axis
end
axis off;

end

