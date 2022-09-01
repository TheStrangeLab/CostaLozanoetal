clear all
close all
clc

restoredefaultpath
addpath 'C:\Users\manuela\Desktop\CTB\000_Toolbox\fieldtrip-20210212\'
ft_defaults

subjects=[5 500 600 8 10];
mlist={'s5','s500','s600','s8','s10'};
behavlist = {'Patient5','Patient500','Patient600','Patient8','Patient10'};

%pz5 left perirhinal
chanofinterest{1} = {'ECL1','ECL2'};
labelorg{1} = {'ECL1','ECL2'};
bipolabel{1} = {'ECL1 - ECL2'};
bipolmat{1} = [+1 -1]; 

% pz5 right perirhinal 
chanofinterest{2} = {'ECR1','ECR2'};
labelorg{2} = {'ECR1','ECR2'};
bipolabel{2} = {'ECR1 - ECR2'};
bipolmat{2} = [+1 -1]; 

% pz6 right perirhinal 
chanofinterest{3} = {'ECR1','ECR2'};
labelorg{3} = {'ECR1','ECR2'};
bipolabel{3} = {'ECR1 - ECR2'};
bipolmat{3} = [+1 -1]; 

% pz8 left perirhinal
chanofinterest{4} = {'ECL1','ECL2'};
labelorg{4} = {'ECL1','ECL2'};
bipolabel{4} = {'ECL1 - ECL2'};
bipolmat{4} = [+1 -1]; 

%pz10 left perirhinal 
chanofinterest{5} = {'ECL1','ECL2'};
labelorg{5} = {'ECL1','ECL2'};
bipolabel{5} = {'ECL1 - ECL2'};
bipolmat{5} = [+1 -1];

%% Paths
oripath = 'C:\Users\manuela\Desktop\AversiveMemFormation';
addpath(fullfile(oripath,'Costalozanoetal','code','utils'))

v=0;
for sub=subjects
 
    %% define subject
v=v+1;

    
    strPaths.Main = ['~\DataZurich\Patient',num2str(sub)];  
%%

load(fullfile(strPaths.Main,'enc_onsets.mat'));

% Main folder of task and subfolders
strPaths.Strange                                  = strPaths.Main; 
% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = '~\fieldtrip-20210212\'
% Matlab Import Import NLX files
strPaths.Toolboxes.MatlabImportExport   = '~\neuralynx-master\MatlabImportExport_v6.0.0\';  
% EEGLAB toolbox (if you have it)
strPaths.Toolboxes.EEGLAB               = '~\eeglab14_1_2b\'; % It shouldn't be on the MATLAB path

% Change main directory
cd(strPaths.Main)

% Add FieldTrip to path
addpath(strPaths.Toolboxes.FieldTrip)
ft_defaults
% Add toolbox for importing NLX to path / This makes importing much faster
addpath(genpath(strPaths.Toolboxes.MatlabImportExport))

% add Get_Clean_Event_Structure to path
addpath(genpath('~\CostaLozanoetal\utils'))

%% Path for the NLX folder
% FieldTrip cannot read the events or header if ncs files do not match in
% their sampling frequency or starting time
% What works is to put macro, scalp and micro data in different folders
% Below is the path for the folder with only 'Events_0003.nev', 'PegasusLogFile.txt' and
% 64 .ncs files for 64 macro channels. NOTHING ELSE.
strNLXFolderPath = ['~\DataZurich\Patient',num2str(sub),'\iaps_enc','\m\']

%% Read header and events
% Don't mind the warnings / "Warning: discontinuous recording, predicted number of timestamps and observed number of timestamps differ by"
hdr = ft_read_header(strNLXFolderPath); 
evR = ft_read_event(strNLXFolderPath);
% CHANGE THIS
% 'Events_0003.nev' (or 'Events_000x.nev', x being any number) should be renamed to 'Events.nev' because that is the file name that FieldTrip searches for
% In a subfunction there is a part for matching the .nev file name to 'Events.nev', 'events.Nev' etc. 
% If the .nev file name doesn't match the templates FieldTrip searches for there will be an error

%% Plot events
ev_ToPlot = evR;
figure
stem(([ev_ToPlot.sample]-1)/hdr.Fs,[ev_ToPlot.value])
ylabel('Trigger value')
xlabel('Time (s)')
xlim([min(([ev_ToPlot.sample]-1))/hdr.Fs,max(([ev_ToPlot.sample]-1))/hdr.Fs])

%% Select the interval for event structure
% Events are read as a struct array, each element corresponding to a trigger
% These are created throughout the recording, not necessarily throughout the experiment
% This includes trial runs, errors etc.
% So the whole length of the event struct array may not be suitable
% Before this is usable, you should determine the range of events to use
% Here I set the range as from the beginning to the end
% Select the range based on the plot of event triggers vs time
indEventRangeStart = 1;%4; % Beginning of event structure % CHANGE THIS
indEventRangeStop = length(evR); % End of event structure % CHANGE THIS

evSelected = evR(indEventRangeStart:indEventRangeStop);

%% Plot selected range for events again
ev_ToPlot = evSelected;
figure
stem(([ev_ToPlot.sample]-1)/hdr.Fs,[ev_ToPlot.value])
ylabel('Trigger value')
xlabel('Time (s)')
xlim([min(([ev_ToPlot.sample]-1))/hdr.Fs,max(([ev_ToPlot.sample]-1))/hdr.Fs])

%% Clean event structure
% The event structure needs to be 'cleaned'
% Normally, every nonzero trigger should be preceded and followed by a zero
% If there are two consecutive nonzero triggers, this is because of the
% rising or falling edge of the trigger signal
% These neeed to be removed by the following function
evClean = Get_Clean_Event_Structure(evSelected);

%% Read timestamps from a single channel
 strFirstFileName = dir(fullfile(strNLXFolderPath,'*.ncs'));
 strFirstFileName = {strFirstFileName.name}';
 strFirstFileName = strFirstFileName(1);
 strFirstFilePath = strcat(strNLXFolderPath,strFirstFileName);

% Timestamp corresponding to each sample
% ts = ft_read_data(strFirstFilePath,'timestamp', true);
ts = ft_read_data(strFirstFilePath,'timestamp', true);
figure,plot(ts)

% Samples correspond to whatever data is saved on disk within an ncs file
% Timestamps, on the other hand, keep increasing independent of the recording
% Normally if there are no gaps in the timing of acquisition and recording,
% samples and timestamps follow the following equation
% timestamp = hdr.FirstTimeStamp + (sample-1)*hdr.TimeStampPerSample;
% The event structure assumes this equation is true
% This will plot the linear relationship
% figure,plot(([evClean.sample]-1)./(double([evClean.timestamp]-hdr.FirstTimeStamp)/hdr.TimeStampPerSample))
% However, when the recording is discontinuous, there are not always
% samples recorded corresponding to timestamps, so the equation does not
% hold

%% Correct samples of the event structure
evCorrected = evClean;
for ii = 1:length(evCorrected)
    % The sample corresponding to this timestamp
    [~,indSample] = min(abs(double(ts)-double(evCorrected(ii).timestamp)));
    % Update field 'sample'
    evCorrected(ii).sample = indSample;
end
% Plot for both event structures
figure,plot(([evClean.sample]-1)./(double([evClean.timestamp]-hdr.FirstTimeStamp)/hdr.TimeStampPerSample))
figure,plot(([evCorrected.sample]-1)./(double([evCorrected.timestamp]-hdr.FirstTimeStamp)/hdr.TimeStampPerSample))

%% Plot all timestamps and timestamps only for the event structure together
% The timestamps should always be monotonically increasing
% If this is not true outside of the time where the experiment is done, no
% problem
% to check if timestamps are monotonically increasing during the experiment,
% I plot here timestamps from the whole recording and from the triggers
% of the experiment together
% I also plot on top the intervals where the timestamp behaves differently
% (green curve)
figure
plot(ts,100*zeros(length(ts),1),'b')
hold on
plot([evClean.timestamp],[evClean.value],'ro','MarkerFaceColor','r')
plot(ts(1:end-1),diff(ts)<0,'g')
legend('Timestamps','Events','Timestamp errors')

%% Define trials 
% Define sample ranges for each trial
% This depends on your task
t_pre_stim = 7.5; % How many seconds before the trigger to take
t_post_stim = 7.5; % How many seconds after the trigger to take

cfg = []; % initialize
cfg.dataset = strNLXFolderPath;
cfg.trialdef.eventtype  = 'trigger'; % Define trials based on triggers
cfg.trialdef.eventvalue = 10; % Select trigger that is repeated only once in each trial
cfg.trialdef.prestim    = t_pre_stim;
cfg.trialdef.poststim   = t_post_stim;
cfg.event = evCorrected; % Event structure for the trials
cfg = ft_definetrial(cfg); % If your triggers are more complex, you might need to write your own functions

% get emotional trials only
    trl = cfg.trl; %remember all trials
    cfg.trl = trl(enc_onsets.eCorrRem,:);

    % now read the data

    cfg.detrend = 'yes';
    cfg.continuous = 'yes';
    cfg.montage.tra      = bipolmat{v};
    cfg.montage.labelnew = bipolabel{v};
    cfg.montage.labelorg = labelorg{v};
    
    data_eR = ft_preprocessing(cfg);
     
    cfg.trl = trl(enc_onsets.eMissed,:);
    data_eF = ft_preprocessing(cfg);
    
    cfg.trl = trl(enc_onsets.nCorrRem,:);
    data_nR = ft_preprocessing(cfg);
    
    cfg.trl = trl(enc_onsets.nCorrFam, :);
    data_nK = ft_preprocessing(cfg);
    
    cfg.trl = trl(enc_onsets.nMissed,:);
    data_nF = ft_preprocessing(cfg);
    
    cfg.trl = trl(enc_onsets.eCorrFam,:);
    data_eK = ft_preprocessing(cfg);
    
   % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eRclean = ft_rejectvisual(cfg,data_eR);
    
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nRclean = ft_rejectvisual(cfg,data_nR);
    
         % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eFclean = ft_rejectvisual(cfg,data_eF);
  
%     % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nFclean = ft_rejectvisual(cfg,data_nF);
% 
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_eKclean = ft_rejectvisual(cfg,data_eK);
     
    % look for spikes
    cfg          = [];
    cfg.channel = bipolabel{v};
    cfg.latency = [-0.5 1.5] 
    cfg.keeptrial = 'nan' %  'nan' fill the channels that are deselected with NaNs
    cfg.method   = 'channel';
    data_nKclean = ft_rejectvisual(cfg,data_nK);
  %% save preprocessed data   
  
% file = [behavlist{v},'_CleanTrialsParaHcECPC_memory_bipolar.mat'];
% save(file,'data_eRclean','data_eKclean','data_eFclean','data_nRclean','data_nKclean','data_nFclean');
  
close all
end
