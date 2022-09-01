function events = runRemoveBlinks(events,outDir,Channels,eogChannels,winsizeMS,overlapMS,corrThresh)
%RUNREMOVEBLINKS - Runs removeBlinks :)
%
% This is a wrapper function for the removeBlinks function, which
% uses a PCA-based algorithm to remove EOG-related artifacts from
% all EEG channels.
%
% You should run this function after detecting bad channels and
% prior to rereferencing the data.
%
% RUNREMOVEBLINKS will not overwrite files, so if you want to rerun
% it, you must first delete the output files it created.
%
% FUNCTION:
%   events = runRemoveBlinks(events,outDir,Channels,eogChannels,winsizeMS,overlapMS,corrThresh)
%
% INPUT ARGS:
%   events
%   outDir
%   Channels
%
%   See removeBlinks function for:
%     eogChannels
%     winsizeMS
%     overlapMS
%     corrThresh
%

% OUTPUT ARGS:
%   events - Events with new eegfile value with outDir.
%
%
% 
%
%

% handle default vals
if ~exist('eogChannels','var')
  eogChannels = [];
end
if ~exist('winsizeMS','var')
  winsizeMS = [];
end
if ~exist('overlapMS','var')
  overlapMS = [];
end
if ~exist('corrThresh','var')
  corrThresh = [];
end

% Get samplerate and data info from event
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(1));

% Loop over unique data files
eegfile = unique(getStructField(events,'eegfile'));
eegfile = eegfile(~strcmp(eegfile,''));

% create outDir if not existing
if ~exist(outDir,'dir')
  [basepath,outDir] = fileparts(outDir);
  mkdir(basepath,outDir);
end

processedFile = 0;
for f = 1:length(eegfile)
  % Tell what file we are on
  fprintf('%s:\n',eegfile{f});
  
  % ignore bad channels for that data file
  badchanfile = fullfile(eegfile{f},'.bad_chan');
  if exist(badchanfile,'file')
    % load in the bad chans
    fid = fopen(badchanfile,'r');
    bad = fscanf(fid,'%d');
    
    % remove from Channels
    Channels = setdiff(Channels,bad);
  end
  
  % see if already in progress
  c = Channels(1);
  [fpath,basename] = fileparts(eegfile{f});
  chanfile = [fullfile(outDir,basename) '.' num2str(c,'%03d')];
  if exist(chanfile,'file')
    % already in progress or done, so continue
    continue;
  end
  
  % touch the file to reserve
  fid = fopen(chanfile,'wb','l');
  fclose(fid);
  processedFile = 1;
    
  % Load in data from all channels
  % make fake event
  event.eegfile = eegfile{f};
  event.eegoffset = 1;
  eeg = [];
  for c = Channels
    fprintf('%d ',c);
    e = gete(c,event,0);
    
    if isempty(eeg)
      eeg = zeros(max(Channels),length(e{1}),'single');
    end
    
    % append the data
    eeg(c,:) = single(e{1});
  end
  fprintf('\n');
    
  % call removeBlinks
  % *** Format is hard-coded as short ***
  outputformat = 'short';
  new_eeg = int16(removeBlinks(eeg,eogChannels,samplerate,winsizeMS,overlapMS,corrThresh)./gain);
  clear eeg;
    
  % write out new data
  for c = Channels
    % set new outfile
    [fpath,basename] = fileparts(eegfile{f});
    chanfile = [fullfile(outDir,basename) '.' num2str(c,'%03d')];
    fid = fopen(chanfile,'wb','l');
    fwrite(fid,new_eeg(c,:),outputformat);
    fclose(fid);
  end
end
  
% write out params.txt file
if processedFile
  paramfile = fullfile(outDir,'params.txt');
  fid = fopen(paramfile,'w');
  fprintf(fid,'samplerate %d\ndataformat ''%s''\ngain %g\n',samplerate,outputformat,gain);
  fclose(fid);
end
           
% % update the eegfile in the events
% for e = 1:length(events)
%   % get the base file
%   [fpath,basename] = fileparts(events(e).eegfile);
%   events(e).eegfile = fullfile(outDir,basename);
% end
  
     
