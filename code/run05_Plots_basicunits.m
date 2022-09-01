%This script reproduce Supplementary Fig.21 a,b,c,d
%% load previously saved results
clear all
close all

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
r   = load(fullfile(oripath,'data2github','SingleUnit','BasicUnitsresults.mat'))
% unit indices
allUnitIdx  = cell2mat({r.allRes.idx}');
%% units per wire

% number of units per wire
uniqueAllWireIdx    = unique(allUnitIdx(:, 1:2), 'rows', 'stable');
fprintf('Number of wires with at least one unit that was analyzed: %d.\n', ...
    size(uniqueAllWireIdx, 1));
numUnitsPerWire     = nan(size(uniqueAllWireIdx, 1), 1);
for iWire = 1:size(uniqueAllWireIdx, 1)
    numUnitsPerWire(iWire)  = sum(all(uniqueAllWireIdx(iWire, 1:2) == allUnitIdx(:, 1:2), 2));
end

% average number of units per wire
fprintf('Number of units per wire: %.3f +/- %.3f (mean +/- SEM).\n', mean(numUnitsPerWire), ...
    std(numUnitsPerWire) / sqrt(size(numUnitsPerWire, 1)));

% histogram showing the number of units per wire
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
histogram(numUnitsPerWire, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'box', 'off', ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel('Units per wire');
yl = ylabel('Number of wires');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
%print(f, strcat(paths.save, 'UnitsPerWire'), '-dtiff', '-r450');

% save data to recreate figure
%save(strcat(paths.save, 'Fig_S3A_data'), 'numUnitsPerWire');

%% ISI refractoriness

% ISI refractoriness from all cells
allPercISILessThan3ms   = cell2mat({r.allRes.percISIlessThan3ms}');

% report
fprintf('Average percentage of ISIs < 3ms across units: %.3f +/- %.3f (mean +/- SEM).\n', ...
    mean(allPercISILessThan3ms), std(allPercISILessThan3ms) / sqrt(size(allPercISILessThan3ms, 1)));
fprintf('Number of units with a percentage of ISIs < 3ms higher than 5%%: %d.\n', ...
    sum(allPercISILessThan3ms >= 5));

% histogram quantifying the ISI refractoryness
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
histogram(allPercISILessThan3ms, 0:0.5:20, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xlim', [0, 5], 'xtick', 0:1:20, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('% ISI <3 ms');
yl = ylabel('Number of units');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
%print(f, strcat(paths.save, 'ISIRefractoryness'), '-dtiff', '-r450');

% save data to recreate figure
%save(strcat(paths.save, 'Fig_S3C_data'), 'allPercISILessThan3ms');

%% firing rate

% mean firing rate from all cells
allMeanFR   = cell2mat({r.allRes.meanFR}');

% report
fprintf('Mean firing rate across all units: %.3f +/- %.3f Hz (mean +/- SEM).\n', ...
    mean(allMeanFR), std(allMeanFR) / sqrt(size(allMeanFR, 1)));

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
histogram(allMeanFR, 0:1:25, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'xtick', 0:5:25, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Mean FR (Hz)');
yl = ylabel('Number of units');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
%print(f, strcat(paths.save, 'MeanFR'), '-dtiff', '-r450');

% save data to recreate figure
%save(strcat(paths.save, 'Fig_S3E_data'), 'allMeanFR');

%% peak SNR

% peak SNR from all cells
allPeakSNR  = cell2mat({r.allRes.peakSNR}');

% report
fprintf('Average peak SNR: %.3f +/- %.3f.\n', mean(allPeakSNR), ...
    std(allPeakSNR) / sqrt(size(allPeakSNR, 1)));

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 8, 8]);
histogram(allPeakSNR, 0:3:33, ...
    'FaceColor', [0.7, 0.7, 0.7]);
set(gca, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Peak SNR');
yl = ylabel('Number of units');
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
set(f, 'PaperPositionMode', 'auto');
%print(f, strcat(paths.save, 'PeakSNR'), '-dtiff', '-r450');

% save data to recreate figure
%save(strcat(paths.save, 'Fig_S3G_data'), 'allPeakSNR');
