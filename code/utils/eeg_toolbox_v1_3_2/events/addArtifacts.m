function events = addArtifacts(eventfile,Channels,thresh,plotit,replace_eegfile)
%ADDARTIFACTS - Add artifact info to events structure.
%
% This function loads an events structure, loops over the unique
% files in the events structure, loads the EEG data from each
% channel in each file, makes the channel data bipolar if
% necessary, finds the blink events for each channel, performs a
% logical OR of the blinks across channels, populates the event
% structure structure with the millisecond time of the first blink
% in each event, then saves the events structure back to file.
%
% In case you need to do more detailed processing of artifact times,
% this function saves the indexes of each artifact to the file
% eegfile_blink.mat for each unique eegfile in the events structure.
% You can load this file with the loadVar function. 
%
% FUNCTION: 
%   events = addArtifacts(eventfile,Channels,thresh,plotit)
%
% INPUT ARGS:
%   eventfile = 'beh/events.mat';  % Events struct to process
%   Channels = {1,2,[63,64]};      % Cell array of channels to
%                                  %  search for eye blinks. Two
%                                  %  element arrays of channels will
%                                  %  have the diff taken before
%                                  %  processing for blinks.
%   thresh = [50,20,50];           % Thresh in uV for each channel
%                                  %  If one value present, will
%                                  %  use for all.
%   plotit = 0;                    % 1 or 0 to plot the signal, signal
%                                  %  distribution, and eyeblinks
%                                  %  for determining a better thresh
%
% OUTPUT ARGS:
%   events - Returns the new events structure that it wrote to file.
%

if ~exist('thresh','var')
  thresh = 100;
end
if ~exist('plotit','var')
  plotit = 0;
end
if ~exist('replace_eegfile','var')
  replace_eegfile = [];
end


% see if we are processing a single file
if ~strcmp(eventfile(end-3:end),'.mat')
  % make fake events
  events = struct('eegfile',{eventfile},'eegoffset',{1});
  
  % get the fnames
  fnames = {eventfile};
  
  % dont process events
  processEvents = 0;
else
  % load the events
  events = loadEvents(eventfile,replace_eegfile);

  % save a backup
  backupfile = [eventfile '.old'];
  saveEvents(events,backupfile);
  fprintf('Saved backup to %s\n',backupfile);

  % get the uniqe file names that are not ''
  fnames = unique(getStructField(events,'eegfile','strcmp(eegfile,'''')==0'));

  % process events
  processEvents = 1;
end

% get data info
tempEvents = filterStruct(events,'strcmp(eegfile,'''')==0');
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(tempEvents(1));

% allocate space for blink index
ind = cell(length(Channels),length(fnames));

% see if using different thresh
if length(thresh) < length(Channels)
  thresh = ones(1,length(Channels)).*thresh(1);
end

% loop over sets of channels
for c = 1:length(Channels)
  % load all the eeg for all the channels
  fprintf('Loading channel data: ');
  for i = 1:length(Channels{c})
    fprintf('%d ',Channels{c}(i));
    eeg{i} = gete(Channels{c}(i),tempEvents,0,0,0,[]);
  end
  fprintf('\n');
  
  % see if must make bipolar
  if length(eeg) > 1
    % must take diff of channels for each file
    for i = 1:length(eeg{1})
      beeg{i} = eeg{1}{i} - eeg{2}{i};
    end
  else
    % just copy the data out of the cell array
    for i = 1:length(eeg{1})
      beeg{i} = eeg{1}{i};
    end
  end
  
  % now get the artifacts for each channel
  fprintf('Searching for blinks using thresh=%guV...\n',thresh(c));
  for i = 1:length(beeg)
    [ind{c,i},fast{i},slow{i}] = findBlinks(beeg{i},thresh(c));
  end
  
  % see if user wanted a plot of it
  if plotit
    for i = 1:length(beeg)
      num_row = length(Channels)*length(fnames);
      cur_plot = ((c-1)*length(fnames))+i;
      cur_plot = ((c-1)*length(fnames)*2)+i;
      subplot(num_row,2,cur_plot)
      %subplot(num_row,1,cur_plot)
      plot(beeg{i})
      pind = find(ind{c,i});
      if length(pind) > 1
	hold on
	plot(pind,max(beeg{i})/2,'r.')
	hold off
      end
      title(['File: ' fnames{i} ', Channels: ' num2str(Channels{c})])
    
      subplot(num_row,2,cur_plot+1);
      % get the max of each window
      %ex_beeg = extend(beeg{i},window);
      %ex_beeg = reshape(ex_beeg,round(length(ex_beeg)/window),window);
      %winmax = squeeze(max(abs(ex_beeg),[],2)-mean(ex_beeg,2));
      %[in,winmax] = local_max(abs(fast{i}));
      %hist(winmax,500);
      drawnow
    end
  end

end

% clear out the eeg data
%clear eeg beeg

% do logical across channels for all possible blinks for each file
blinks = cell(1,size(ind,2));
for i = 1:length(blinks)
  blinks{i} = ind{1,i};
  for c = 2:length(Channels)
    blinks{i} = or(blinks{i},ind{c,i});
  end

  % get blink indexes
  blinks{i} = find(blinks{i});
end

% see if we have events to process
if ~processEvents
  % return the blinks and do no more
  events = blinks{1};
  return
end

% write out blink indexes
for i = 1:length(blinks)
  blinkfile = [fnames{i} '_blinks.mat'];
  saveVar(blinks{i},blinkfile);
end

% now we have the blink indexes, loop over events
fprintf('Adding artifacts to events...\n');
for e = 1:length(events)
  % make sure the event has an eegfile
  if length(events(e).eegfile) == 0
    % no file, so skip it adding the field with no blink
    events(e).artifactMS = -1;
    continue
  end
  
  % see which file we are in
  i = find(ismember(events(e).eegfile,fnames));
  
  % look for blink
  if e+1 < length(events)
    % can use next event as upper bound
    ev_blink = blinks{i}(blinks{i}>=events(e).eegoffset & blinks{i} < events(e+1).eegoffset);
  else
    % no upper bound
    ev_blink = blinks{i}(blinks{i}>=events(e).eegoffset);
  end
  
  % if have blink, convert first one to ms
  if ~isempty(ev_blink)
    % add in the first one
    blinkOffset = (ev_blink(1) - events(e).eegoffset);
    blinkMS = blinkOffset*1000/samplerate;
    events(e).artifactMS = blinkMS;
  else
    % no blink
    events(e).artifactMS = -1;
  end
  
end

% save the updated events
saveEvents(events,eventfile);
fprintf('Saved updated events to %s.\n',eventfile);



