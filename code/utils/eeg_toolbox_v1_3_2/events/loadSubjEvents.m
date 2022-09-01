function events = loadSubjEvents(basedir,subj,eventfile,replace_eegfile)
%LOADSUBJEVENTS - Combine a number of subject's events
%
% Given a base directory and list of subjects, this function will
% loop through and load each subject's events and concatenate them 
% into a single events structure.
%
% FUNCTION:
%   events = loadSubjEvents(basedir,subj,eventfile,replace_eegfile)
%
% INPUT ARGS:
%   basedir = '~/eeg/free';   % Root directory to look for files
%   subj = {'CH003','CH005'}; % Cell array of subjects to load
%   eventfile = 'events/events.mat';  % directory and file to load
%   replace_eegfile = {'/Users/lynne/FRdata','/data/eeg/scalp/fr/fr1'};
%
% OUTPUT ARGS:
%   events- Events struture of all combined events
%

% 05/29/05 - pbs - replace_eegfile now can prepend if you give just one.
% 05/10/05 - pbs - Now does not go into directory to load.
% 03/23/05 - pbs - Added replace_eegfile to enable changing of eegfile
%                  on load.
% 10/29/04 - pbs - Made it so it adds the fields to match the structures.


events = [];

% make sure we have cells
if ~iscell(subj)
  subj = {subj};
end
if ~iscell(eventfile)
  eventfile = {eventfile};
end

% loop over subjects
for s = 1:length(subj)
  % go to 
  subjdir = fullfile(basedir,subj{s});
  
  % loop over eventfiles
  for e = 1:length(eventfile)
    % create file
    efile = fullfile(subjdir,eventfile{e});
    
    % make sure file exists
    if ~exist(efile,'file')
      continue
    end
    
    % see if must pick a replace_eegfile
    if exist('replace_eegfile','var') & ~isempty(replace_eegfile)
      % must handle eegfile
      if size(replace_eegfile,1) > 1
        % pick one associated with current eventfile
        current_replace_eegfile = {replace_eegfile{e,:}};
      else
        % just pick first one
        current_replace_eegfile = replace_eegfile;
      end
      
    else
      current_replace_eegfile = [];
    end
    
    
    % load the events
    new_events = loadEvents(efile,current_replace_eegfile);
    
    if ~isempty(events)
      % add dummy fieldnames for those that don't match
      addtonew = setdiff(fieldnames(events),fieldnames(new_events));
      for i = 1:length(addtonew)
        new_events(1).(deblank(addtonew{i})) = [];
      end
      addtoold = setdiff(fieldnames(new_events),fieldnames(events));
      for i = 1:length(addtoold)
        events(1).(deblank(addtoold{i})) = [];
      end
    end
    
    
    % concatenate the structures
    events = [events , new_events];
  end
  
end

