function events = createEvents(event_file,event_format,event_fields,varargin)
%CREATEEVENTS - Creates an events structure from a file
%
% Create an events structure from a generic text file.
% You can then use the functions getStructField and filterStruct to query for
% different combinations of events.
%
% The basis of the EEG_TOOLBOX is the event structure.  Usually,
% the person who designs an experimental paradigm will also write
% a custom createEvents function that will create the events
% structure from the experiment-specific log file.  In cases where a user
% would like to simply import event information so that they can
% use the toolbox, we provide this generic function that will
% create the events structure needed to use most of the analysis
% functions.
%
% The main requirement for the EEG_TOOLBOX events structure is that
% it have an 'eegfile' field and an 'eegoffset' field.  The
% 'eegfile' field defines the path and root of the eegfiles:
% '/data/eeg/CH012/dat/CH012_13May03_1203'.  The analysis functions
% will automatically add the channel suffix.  The 'eegoffset' field
% defines the offset in samples into the eeg data file when each
% event occured.  The functions that align your eeg and behavioral
% data will automatically add the eegfile and eegoffset fields
% based on your eeg sync pulses and the mstime field.
%
%
% FUNCTION: 
%   events = createEvents(event_file,event_format,event_fields)
%
%
% INPUT ARGS:
%   event_file = 'session.log';  % Text event file to read in
%   event_format = '%n%n%s%n';  % Column format in event file
%   event_fields = {'mstime','msoffset','item','recalled'}; % Name of the fields
%   varargin;     % Any additional args that are passed directly to textread.
%
% OUTPUT ARGS:
%   events- An events structure that you should save to
%           events/events.mat with saveEvents and load when needed
%           with loadEvents.
%


if ~exist('varargin','var')
  varargin = {};
end

% create string for loading the columns
colstr = '[';
for i = 1:length(event_fields)
  colstr = [colstr 'field_' num2str(i)];
  if i < length(event_fields)
    colstr = [colstr ','];
  end
end
colstr = [colstr ']'];

% load fields from the file
strtext = [colstr ' = textread(event_file,event_format,varargin{:});'];
eval(strtext);

% now make the struct
strstruct = ['events = struct('];

for i = 1:length(event_fields)
  eval(['addcell = ~iscell(field_' num2str(i) ');']);
  if addcell
    eval(['field_' num2str(i) ' = num2cell(field_' num2str(i) ');']);
  end
  
  if i>1
    % add the comma
    strstruct = [strstruct ','];
  end
  
  strstruct = [strstruct '''' event_fields{i} ''',field_' num2str(i)];    
end
strstruct = [strstruct ');'];

eval(strstruct);

