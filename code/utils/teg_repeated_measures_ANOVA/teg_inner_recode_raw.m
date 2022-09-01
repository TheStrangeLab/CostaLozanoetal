function [y_red_M, expanded0] = teg_inner_recode_raw(M_raw, NM, cellsets, iPred)

y_red_M = [];
expanded0 = zeros(size(M_raw));
for iCellSet = 1:length(cellsets{iPred}),
    y_raw_cell_obs = M_raw(:, cellsets{iPred}{iCellSet});
    % Take differing cell counts into account
    y_cell_obsN = NM(:, cellsets{iPred}{iCellSet});
    y_cell_obsMN = sum(y_cell_obsN, 2); % Get mean per subject over variable combinations.
    y_cell_obsProp = y_cell_obsN ./ (y_cell_obsMN * ones(1, size(y_cell_obsN, 2)));
    y_raw_cell_obs = y_raw_cell_obs .* y_cell_obsProp;
    tmp = [];
    for iSubj = 1:size(y_raw_cell_obs, 1),
        fNotEmpty = ~isnan(y_raw_cell_obs(iSubj, :));
        tmp(iSubj, 1) = sum(y_raw_cell_obs(iSubj, fNotEmpty), 2);
    end;
    y_raw_cell_obs = tmp;
    y_red_M = [y_red_M y_raw_cell_obs];
    expanded0(:, cellsets{iPred}{iCellSet}) = y_raw_cell_obs * ones(1, length(cellsets{iPred}{iCellSet}));
end;
