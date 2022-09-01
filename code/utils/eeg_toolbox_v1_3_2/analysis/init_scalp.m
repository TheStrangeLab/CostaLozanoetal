function eeg = init_scalp(dataroot, resDir, sessions, params, experiment, anaType, elecLocsFile)
%
%INIT_SCALP - after post-processing, running this script prepares
%one scalp subject's data for analysis
%
% FUNCTION: eeg = init_scalp(dataroot, resDir, sessions, params, experiment, anaType)
%
% INPUT: dataroot - directory containing this subject's data
%        resDir - directory in which to save eeg results
%        sessions - filename of m-file that outputs a 'subj'
%                   struct, or the subj struct itself (see README).
%        params - struct containing eeg analysis parameters (see
%                 below for default settings, and README for
%                 possible additional fields)
%        experiment - name of the experiment (optional)
%        anaType - name of the analysis to be carried out (optional)
%
% OUTPUT: eeg, a struct containing all basic info for this subject;
% gets passed into all other eeg analysis scripts
%
% EXAMPLE:
% dataroot = '/data/eeg/scalp/catFR/subj_00';
% resDir = '/users/morton/EXPERIMENTS/catFR/results';
% sessions = 'iCatFR_sessions.m';
%            or a subj struct: subj(1).id = subj_00; subj.sess(1).eventsFile = '/data/eeg/scalp/catFR/subj_00/session_0/events.mat'
% params = struct('durationMS', 1800, 'offsetMS', -200...);
% experiment = 'catFR';
% anaType = 'power';

% set minimum defaults
if nargin<7
  elecLocsFile = '';
  if nargin<6
    anaType = 'unknown';
    if nargin<5
      experiment = 'unknown';
      if nargin<4
	params = struct('channels', 1:129);
      end
    end
  end
end

if ~exist(resDir)
  mkdir(resDir);
end  

% create the eeg struct
eeg = struct('experiment', experiment, 'recordingType', 'scalp', 'anaType', anaType, 'dataroot', dataroot, 'resDir', resDir, 'params', params);

% add eventsFile info for each subj, session
if isstr(sessions)
  run(sessions);
  eeg.subj = subj;
elseif isstruct(sessions)
  eeg.subj = sessions;
end

% if an electrode locations file is available, read it
if ~isempty(elecLocsFile) & exist(elecLocsFile, 'file')
  [channels regions] = textread(elecLocsFile, '%d%s'); 
end

for c=1:length(params.channels)
  eeg.chan(c).number = params.channels(c);
  if exist('regions','var')
    eeg.chan(c).region = regions{channels==eeg.chan(c).number};
  else
    eeg.chan(c).region = '';
  end
end

for s=1:length(eeg.subj)
  
  % each subject gets the same channel info
  eeg.subj(s).chan = eeg.chan;
  
  % for each session, find out which channels were good
  for n=1:length(eeg.subj(s).sess)
    bad_chan_dir = fullfile(fileparts(eeg.subj(s).sess(n).eventsFile), 'eeg');
    temp = dir(fullfile(bad_chan_dir, '*.bad_chan'));
    bad_chans = textread(fullfile(bad_chan_dir, temp.name));
    eeg.subj(s).sess(n).goodChans = setdiff(params.channels, bad_chans);
  end
end

save(fullfile(eeg.resDir, 'eeg.mat'), 'eeg');


