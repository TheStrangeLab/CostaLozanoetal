%==========================================================================
% This script performs spike-detection and -clustering using wave-clus.
% Lukas Kunz, 2021/01/30
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('~\wave_clus-NEW')); % use new wave-clus (Chaure et al., 2018)

% spike direction to detect
param           = [];
param.spikeDir  = 'neg'; % neg | pos | both

% paths
paths           = [];
paths.data      = '~\DataZurich\Micro\';
paths.spike     = strcat('C:\Users\manuela\Desktop\AversiveMemFormation\data2github\SingleUnit\SpikeExtraction_20210210\', param.spikeDir, 'Spikes\');
paths.pic       = strcat(paths.spike, 'AllPics\');
if ~exist(paths.pic, 'dir')
    mkdir(paths.pic);
end

% subjects
subjects    = {...
    'Patient1'; ...
    'Patient2'; ...
    'Patient5'; ...
    'Patient6'; ...
    'Patient8'; ...
    'Patient10'; ...
    };

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report progress
    rng(1); % for reproducibility
    fprintf('\nWorking on subject %s ...\n', subjects{iSub});
    cd(paths.spike);
        
    % wires to process
    files   = dir(strcat(paths.data, '\', subjects{iSub}, '\u*'));    
    fprintf('\n------------ Performing wave_clus\n');
        
    %% loop through wires and extract action potentials
    for iFile = 1:size(files, 1)
            
        fprintf('\tWorking on wire "%s" ...\n', fullfile(files(iFile).folder, files(iFile).name));
            
        % storage destination
        paths.spike_spec = fullfile(paths.spike, subjects{iSub}, filesep, files(iFile).name, filesep);
            
        % copy original data to storage destination
        mkdir(paths.spike_spec);
        status = copyfile(fullfile(files(iFile).folder, files(iFile).name, filesep, 'data.mat'), ...
            paths.spike_spec);
            
        %% wave clus
        % needs input file containing "data" and "sr"
        
        % get file for wave-clus
        cd(paths.spike_spec);
        file2use4wc     = {strcat(paths.spike_spec, 'data.mat')};
        
        % run wave_clus to extract spikes and do clustering
        try
            % perform spike extraction
            myPar                   = [];
            myPar.detection         = param.spikeDir; % determine detection type
            myPar.randomseed        = 1; % define randomseed
            Get_spikes(file2use4wc, 'par', myPar);
            
            % perform spike clustering
            file2use4wc             = {strcat(paths.spike_spec, 'data_spikes.mat')};
            myPar                   = [];
            myPar.randomseed        = 1; % define randomseed
            myPar.template_sdnum    = 1.5; % for each unsorted spike, look up the closest template - but only if it is within its SD times 'par.template_sdnum'; (default 3)
            myPar.min_clus          = 60; % minimum size of a cluster; (default 20)
            myPar.max_clus        	= 10; % maximum number of clusters allowed (default 13)
            myPar.mintemp           = 0.05; % minimum temperature
            Do_clustering(file2use4wc, 'par', myPar);
            
            % copy output picture(s) to overview folder
            pics    = dir(fullfile(paths.spike_spec, '*.png'));
            for iPic = 1:length(pics)
                copyfile(fullfile(pics(iPic).folder, pics(iPic).name), ...
                    strcat(paths.pic, subjects{iSub}, '_', files(iFile).name, '_', pics(iPic).name));
            end
        catch
            warning('Wave-clus could not be performed.');
        end
    end
end