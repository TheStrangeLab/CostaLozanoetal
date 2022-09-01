%==========================================================================
% This script extracts the microwire data from the original files for use
% with wave-clus.
% Lukas Kunz, 2021/01/30
%==========================================================================

% start
clc; close all; clear;

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20201113\'
ft_defaults

oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))
addpath(fullfile(oripath,'Costalozanoetal','code','utils','eeg_toolbox_v1_3_2'))% kahana eeg toolbox
addpath(fullfile(oripath,'Costalozanoetal','code','utils','MatlabImportExport_v6.0.0'));
warning('off', 'all');% suppress all warnings


% paths
paths           = [];
paths.micro     = '~\DataZurich\';
paths.save      = '~\DataZurich\Micro\';

% subjects and channels to process
subjects    = {...
    'Patient1', {'uAL'; 'uAHL'; 'uAR'; 'uAHR'}; ...
    'Patient2', {'uAL'; 'uHL'; 'uAR'; 'uHR'}; ...
    'Patient5', {'uAL'; 'uAHL'; 'uAR'; 'uAHR'}; ...
    'Patient6', {'uAL'; 'uHL'; 'uAR'; 'uAHR'}; ...
    'Patient8', {'uAL'; 'uAHL'; 'uAR'; 'uAHR'}; ...
    'Patient10', {'uAL'; 'uAHL'; 'uAR'; 'uAHR'}; ...
    };

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    %% loop through regions
    for iReg = 1:size(subjects{iSub, 2}, 1)
        
        % available files for this patient and region
        files   = dir(strcat(paths.micro, subjects{iSub, 1}, '\EEGData\', subjects{iSub, 2}{iReg}, '*'));
        
        %% loop through files and read them
        for iFile = 1:size(files, 1)
            
            % read neuralynx data using fieldtrip
            cfg         = [];
            cfg.dataset = fullfile(files(iFile).folder, files(iFile).name);
            ncs         = ft_preprocessing(cfg);
            
            % get relevant information for waveclus
            sr          = ncs.fsample;
            data        = ncs.trial{1};
            
            % save data in waveclus format
            tmp             = split(files(iFile).name, '.');
            chanName        = tmp{1};
            thisSavePath    = strcat(paths.save, subjects{iSub, 1}, '\', chanName);
            mkdir(thisSavePath);
            save(strcat(thisSavePath, '\data.mat'), 'data', 'sr');            
        end        
    end
end