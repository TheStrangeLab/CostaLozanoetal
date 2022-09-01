function eeg = create_volt_pattern(eeg, resDir)
% eeg = create_volt_pattern(eeg, resDir)
%
% create a voltage pattern for each subject, time bin
% 


% set the defaults for params
eeg.params = structDefaults(eeg.params,  'eventFilter', '',  'offsetMS', -200,  'durationMS', 1800,  'binSizeMS', 50,  'baseEventFilter', '',  'baseOffsetMS', -200,  'baseDurationMS', 200,  'filttype', 'stop',  'filtfreq', [58 62],  'filtorder', 4,  'bufferMS', 1000,  'resampledRate', 500,  'kthresh', 5,  'ztransform', 1,  'replace_eegFile', {},  'mask', []);

% get bin information
durationSamp = fix(eeg.params.durationMS*eeg.params.resampledRate./1000);
binSizeSamp = fix(eeg.params.binSizeMS*eeg.params.resampledRate./1000);
nBins = fix(durationSamp/binSizeSamp);

binSamp{1} = [1:binSizeSamp];
for b=2:nBins
  binSamp{b} = binSamp{b-1} + binSizeSamp;
end

% get MS values
for b=1:length(binSamp)
  eeg.params.binMS{b} = fix((binSamp{b}-1)*1000/eeg.params.resampledRate) + eeg.params.offsetMS;
end

params = eeg.params;
disp(params);

% prepare dir for the patterns
if ~exist(fullfile(eeg.resDir, 'data'), 'dir')
  mkdir(fullfile(eeg.resDir, 'data'));
end

for s=1:length(eeg.subj)
  
  % see if this subject has been done
  eeg.subj(s).patFile = fullfile(eeg.resDir, 'data', [eeg.subj(s).id '_voltpat.mat']);
  if ~lockFile(eeg.subj(s).patFile)
    save(fullfile(eeg.resDir, 'eeg.mat'), 'eeg');
    continue
  end

  % get all events for this subject, w/filter that will be used to get voltage
  nEvents = 0;
  for n=1:length(eeg.subj(s).sess)
    temp = loadEvents(eeg.subj(s).sess(n).eventsFile, params.replace_eegFile);
    events{n} = filterStruct(temp, params.eventFilter);
    base_events{n} = filterStruct(temp, params.baseEventFilter);
    nEvents = nEvents + length(events{n});
  end 
  
  % check if using custom channels
  if isfield(params, 'channels')
    channels = params.channels;
  else
    channels = getStructField(eeg.subj(s).chan, 'number');
  end
  
  % initialize this subject's pattern
  patSize = [nEvents, length(channels), length(params.binMS)];
  pat.name = eeg.subj(s).id;
  pat.mat = NaN(patSize);
  
  % set up masks
  m = 1;
  mask(m).name = 'kurtosis';
  mask(m).mat = false(patSize);
  m = m + 1;
  mask(m).name = 'bad_channels';
  mask(m).mat = false(patSize);
  if isfield(params, 'artWindow')
    m = m + 1;
    mask(m).name = 'artifacts';
    mask(m).mat = false(patSize);
  end
  for i=1:length(params.mask)
    mask(m+i).name = params.mask(i).name;
    mask(m+i).mat = false(patSize);
  end
  
  % make masks
  e = 1;
  for n=1:length(eeg.subj(s).sess)
    
    % bad channels
    bad = setdiff(channels, eeg.subj(s).sess(n).goodChans);
    for sess_e=1:length(events{n})
      mask(2).mat(e,bad,:) = 1;
      
      if isfield(params, 'artWindow')
	% blink artifacts
	wind = [events{n}(sess_e).artifactMS events{n}(sess_e).artifactMS+params.artWindow];
	isArt = 0;
	for b=1:length(params.binMS)
	  if wind(1)>params.binMS{b}(1) & wind(1)<params.binMS{b}(end)
	    isArt = 1;
	  end
	  mask(3).mat(e,:,b) = isArt;
	  if isArt & wind(2)<params.binMS{b}(end)
	    isArt = 0;
	  end
	end
      end
      
      % custom masks
      for i=1:length(params.mask)
	mask(m+i).mat(e,:,:) = inStruct(events{n}(sess_e), params.mask(i).filter);
      end
      
      e = e + 1;    
    end
  end
  
  start_e = 1;
  for n=1:length(eeg.subj(s).sess)
    fprintf('\n%s\n', eeg.subj(s).sess(n).eventsFile);
    
    for c=1:length(channels)
      fprintf('%d.', channels(c));
      
      % get baseline stats for this channel, sess
      if params.ztransform
	base_eeg = gete_ms(channels(c), base_events{n}, ...
	                   params.baseDurationMS, ...
			   params.baseOffsetMS, params.bufferMS, ... 
			   params.filtfreq, params.filttype, ...
			   params.filtorder, params.resampledRate, ...
			   params.baseRelativeMS);
	
	if ~isempty(params.kthresh)
	  base_eeg = run_kurtosis(base_eeg, params.kthresh);
	end
	  
	poss = 1:size(base_eeg,2);
	for e=1:size(base_eeg,1)
	  randposs = randperm(length(poss));
	  rand_eeg(e) = base_eeg(e,randposs(1));
	end
	
	base_mean = nanmean(rand_eeg(:));
	base_std = nanstd(rand_eeg(:));	   
      end
      
      % get power, z-transform, average each time bin
      e = start_e;
      for sess_e=1:length(events{n})
	this_eeg = squeeze(gete_ms(channels(c), events{n}(sess_e), ...
	                   params.durationMS, params.offsetMS, ...
			   params.bufferMS, params.filtfreq, ...
			   params.filttype, params.filtorder, ...
			   params.resampledRate, params.relativeMS));
	
	if ~isempty(params.kthresh)
	  kInd = kurtosis(this_eeg)>params.kthresh;
	else 
	  kInd = 0;
	end
	mask(1).mat(e,c,:) = kInd;
	
	if params.ztransform
	  this_eeg = (this_eeg - base_mean)/base_std;
	end
	
	if ~isempty(this_eeg)
	  for b=1:nBins
	    pat.mat(e,c,b) = nanmean(this_eeg(binSamp{b}));
	  end
	end
	
	e = e + 1;
      end % events
      
    end % channel
    start_e = start_e + length(events{n});
    
  end % session
  fprintf('\n');
    
  % save the patterns, selectors, and regressors
  save(eeg.subj(s).patFile, 'pat', 'mask');
  releaseFile(eeg.subj(s).patFile);
  save(fullfile(eeg.resDir, 'eeg.mat'), 'eeg');

end % subj




