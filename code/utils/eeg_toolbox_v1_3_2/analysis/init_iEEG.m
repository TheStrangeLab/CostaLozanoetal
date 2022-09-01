function eeg = init_iEEG(dataroot, resDir, sessions, params, experiment, anaType)
%
%INIT_IEEG - after post-processing, running this script prepares
%one iEEG subject's data for analysis
%
% FUNCTION: eeg = init_iEEG(dataroot, resDir, sessions, params, experiment, anaType)
%
% INPUT: dataroot - directory containing this subject's data
%        resDir - directory in which to save eeg results
%        sessions - filename of m-file that outputs a 'subj' struct:
%                       subj.id
%                           .sess(n).eventFile
%                   (or you can pass in the struct itself)
%        params - struct containing eeg analysis parameters (see
%                 below for default settings). possible
%                 additional fields: resampledRate = 500 artThresh =
%                 1600 kthresh = 5 width = 6
%        experiment - name of the experiment (optional)
%        anaType - name of the analysis to be carried out (optional)
%
% OUTPUT: eeg, a struct containing all basic info for this subject;
% gets passed into all other eeg analysis scripts
%
% EXAMPLE:
% dataroot = '/data/eeg/UP011';
% resDir = '/users/morton/EXPERIMENTS/iCatFR/results';
% sessions = 'iCatFR_sessions.m';
%            or subj.id = UP011; subj.sess(1).eventsFile = '/data/eeg/UP011/behavioral/catFR/session_0/events.mat'
% params = struct('durationMS', 1800, 'offsetMS', -200...);
% experiment = 'iCatFR';
% anaType = 'power';
%

% set minimum defaults
if nargin<6
  anaType = 'unknown';
  if nargin<5
    experiment = 'unknown';
    if nargin<4
      params = struct('eventFilter', '');
    end
  end
end

% create the eeg struct
eeg = struct('experiment', experiment, 'recordingType', 'iEEG', 'anaType', anaType, 'dataroot', dataroot, 'resDir', resDir, 'params', params);

% add eventsFile info for each subj, session
if isstr(sessions)
  run(sessions);
  eeg.subj = subj;
elseif isstruct(sessions)
  eeg.subj = sessions;
end

for s=1:length(eeg.subj)
  
  % find out which ones were good, find out what brain region each was in
  good_chans_file = fullfile(dataroot, eeg.subj(s).id, 'tal', 'good_leads.txt');
  good_chans = textread(good_chans_file, '%n');
  jacksheet = fullfile(dataroot, eeg.subj(s).id, 'docs', 'jacksheet.txt');
  [channels, regions] = textread(jacksheet, '%d%s');
  
  [channels, gidx, cidx] = intersect(good_chans, channels);
  for c=1:length(channels)
    eeg.subj(s).chan(c).number = channels(c);
    eeg.subj(s).chan(c).region = regions{cidx(c)};
  end
  
  for n=1:length(eeg.subj(s).sess)
    eeg.subj(s).sess(n).goodChans = good_chans;
  end
  
end

% save the struct, which holds filenames of all saved data
if ~exist(eeg.resDir)
  mkdir(eeg.resDir);
end   
save(fullfile(resDir, 'eeg.mat'), 'eeg');

