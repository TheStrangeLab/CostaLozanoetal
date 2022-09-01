function prep_egi_data(subject,session,eventfiles,badchan,ms_field)
%PREP_EGI_DATA - Prepare data collected with EGI system and pyepl.
%
% FUNCTION:
%   prep_egi_data(subject,session,eventfiles,badchan,ms_field)
%
% INPUT ARGS:
%   subject = 'FR062';
%   session = '.';  % set that when called from session dir
%   eventfiles = {'events.mat','revents.mat'};
%   badchan = [1 4 8];
%   ms_field = 'mstime';
%
% This function expects that the events structures are in the
% session_? directory and the raw eeg data are in the session_?/eeg
% directory.  

if ~exist('ms_field','var')
  ms_field = 'mstime';
end

if ~exist('badchan','var')
  badchan = [];
end

if isstr(session)
  sessdir = session;
else
  sessdir = fullfile('../data',subject,['session_' num2str(session-1)]);
end

eegdir = fullfile(sessdir,'eeg');
norerefdir = fullfile(eegdir,'eeg.noreref');
rerefdir = fullfile(eegdir,'eeg.reref');

% process the eeg.eeglog file
fixEEGLog(fullfile(sessdir,'eeg.eeglog'),fullfile(sessdir,'eeg.eeglog.up'));

% extract the raw data
% find file to extract
rawfile = dir(fullfile(eegdir,'*.raw'));
if length(rawfile) == 0
  error('No raw eeg file found.');
  return
end
basename = egi_split(fullfile(eegdir,rawfile.name),subject,norerefdir);

% get the rate info
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(norerefdir);

% align data
runAlign(samplerate,...
         {'eeg.eeglog.up'},...
         {fullfile(norerefdir,[basename '.DIN'])},...
         {fullfile(norerefdir,[basename '.001'])},...
         eventfiles,...
         ms_field,...
         0,1);

% add artifacts to each event file
for i = 1:length(eventfiles)
  addArtifacts(eventfiles{i},{ [26,127] , [8,126]}, 100,0);
end

% deal with bad channels
% see if badchan is a file to open, or to save
if ~isempty(badchan)
  if isstr(badchan)
    % open from file
    badchanfile = badchan;
    badchan = textread(badchanfile)';
  else
    % save to file
    badchanfile = fullfile(eegdir,[basename '.bad_chan']);
    chan2file(badchanfile,badchan);
  end
end
% set to bad plus perif chans
badchan = unique([badchan 127 126 17 128 125 120 44 49 56 63 69 74 82 89 95 100 108 114]);

% rereference the data
% read in and combine all events structures
fileroots = {};
channels = 1:129;
for i = 1:length(eventfiles)
  ev = loadEvents(eventfiles{i});
  newroots = unique(getStructField(ev,'eegfile','~strcmp(eegfile,'''')'));
  fileroots(length(fileroots)+1:length(fileroots)+length(newroots)) = newroots;
end
fileroots = unique(fileroots);
reref(fileroots,channels,rerefdir,{1:129,setdiff(channels,badchan)});

% redo the events to point to the rereferenced data
for i = 1:length(eventfiles)
  % save backup of old eventfile
  copyfile(eventfiles{i},[eventfiles{i} '.noreref'],'f');
  
  % load in the new file, replacing the eegfile
  ev = loadEvents(eventfiles{i},{rerefdir});
  
  % save out again
  saveEvents(ev,eventfiles{i});
end




function fixEEGLog(infile,outfile)
%FIXEEGLOG - Fix pyepl EEG Logfile leaving only UP pulses.
%
%
% FUNCTION:
%   fixEEGLog(infile,outfile)
%
% INPUT ARGS:
%   infile = 'eeg.eeglog';
%   outfile = 'eeg.eeglog.up';
%

% read in the logfile
[mstime, maxoffset, type] = textread(infile,'%s%n%s%*[^\n]','emptyvalue',NaN);

% write out new file
fid = fopen(outfile,'w');
for i = 1:length(type)
  if strcmp(type{i},'UP')
    % save it to file
    fprintf(fid,'%s\t%d\t%s\n',mstime{i},maxoffset(i),type{i});
  end
end
fclose(fid);


function chan2file(filename,chans)
%CHAN2FILE - Write channels to a file
%
%
% FUNCTION:
%   chan2file(filename,chans)
%
% INPUT ARGS:
%   filename
%   chans
%
%


fid = fopen(filename,'w');
fprintf(fid,'%d\n',chans);
fclose(fid);


