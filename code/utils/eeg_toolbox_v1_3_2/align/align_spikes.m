function events = align_spikes(ncsDir,eventsFile,beh_file,outDir)
% FUNCTION events = align_spikes(ncsDir,events,beh_file)
% this function aligns the spike data that has been created using
% semi-automatic cluster cutting with wave_clus (see CML wiki).
% INPUT Args:
% ncsDir = '/data/eeg/UP011/raw/UP011_18Jul07_1149neuralynx' -
% directory with the raw and cluster cut Neuralynx data
% eventsFile = '/data/eeg/UP011/beh/pyFGS/events.mat - file with
% the events structure
% beh_file =
% '/data/eeg/UP011/behavioral/pyFGS/session_0/eeg.eeglog' - file(s)
% with the records of when the sync pulses were sent by the testing
% computer
% outDir = directory in which the cell vectors will be stored
% OUTPUT Args:
% events with a new field: ncs_offset, as well as a directory with
% one binary file for every cell, indicating for each sample
% whether there is a spike at that time
% 
% NOTE THAT AT THIS TIME, MULTIPLE SESSIONS FOR ONE PATIENT ARE NOT
% YET SUPPORTED@!

events = loadEvents(eventsFile);
if ~exist('outDir','var')
  outDir = 'spikes';
end

maxChan = 80; % the maximum channel number (can this differ per
              % patient)

% 1) load neuralynx sync pulses
cd(ncsDir)
if ~exist('Events.nev','file')
  disp('no nev file found\n');
else
  [upstrokes,downstrokes]=load_nev('Events.nev');
  % the events are times in microseconds since the neuralynx
  % machine was turned on
end
% turn the sync pulses in milliseconds (from microseconds)
upstrokes_ms = upstrokes/1000;
% 2) load behavioral sync pulses, which are in ms
beh_ms = textread(beh_file,'%n%*[^\n]','delimiter','\t');
% 3) perform the regression of neuralynx time on behavioral syncs
% (analogous to the code in runAlign.m)
threshMS = 10; mywin = 100;
% first find a subset of pulses that actually correspond by moving
% a window back and forth 
[beh_in,ncs_in]=pulsealign(beh_ms,upstrokes_ms,[],threshMS,mywin,0,1);
% beh_in and ncs_in are matched sync pulses, and can be directly
% regressed to each other to obtain a slope and offset

% find the slope (should be about 1) and offset
% because of limited precision of regression, we have to divide
% both numbers by 10^6
[b,bint,r,rin,stats]=regress((ncs_in/(10^6))',[ones(length(beh_in),1) beh_in/(10^6)]);
% to get the actual value for the intercept, we multiple again by
% 10^6
b(1) = b(1)*(10^6);
% compute maximum deviation
pred = beh_in*b(2)+b(1);
maxdev = max(abs(pred-ncs_in'));
% report statistics
fprintf('%s:\n',ncsDir);
fprintf('\tMax. Dev. = %f\n',maxdev);
fprintf('\tR^2  = %f\n',stats(1));
fprintf('\tSlope = %f\n',b(2));
fprintf('\n');

% 4) loop through all channels, and for those channels that have
% units, make a binary file with 1 when there is a spike at that
% sample, and 0 if not

%create the directory for writing the binary files with spike data
subject = events(1).subject;
basedir = strtok(ncsDir,subject);
completeOutdir = fullfile(basedir,subject,outDir);
if ~exist(completeOutdir,'dir')
  mkdir(completeOutdir);
end
currentCell = 0;
for c = 1:maxChan
  fprintf('%d ',c);
  filename=fullfile(ncsDir,sprintf('times_CSC%d.mat',c));
  disp(['loading channel ' num2str(c) '. Clusters:']);
  
  if ~exist(filename,'file')
    continue;
  end
  load(filename);
  if isempty(cluster_class)
    disp('no spikes found!');
    continue 
  end
  
  if size(cluster_class)==1
    disp('empty cluster_class variable.')
    continue
  end
  
  % it is possible that the ncsFile is bzipped, hence check that)
  ncsFile = fullfile(ncsDir,['CSC' num2str(c) '.ncs'] );
  if exist(ncsFile,'file')
    ncsStats=stat_ncs(ncsFile);
  elseif exist([ncsFile '.bz2'],'file')
    eval(['!bunzip2 ' ncsFile '.bz2']);
    ncsStats=stat_ncs(ncsFile);
    eval(['!bzip2 ' ncsFile]);
  else
    disp('cannot find NCS file\n');
  end
  [x,outFileStem] = fileparts(ncsDir);
  ind = findstr(outFileStem,'neuralynx');
  outFileStem = outFileStem(1:(ind-1));
  dataLength = fix((ncsStats.lastSampleTime-ncsStats.firstSampleTime)/ncsStats.sampleTSDiff);
  if c==1
    samplerate = 1/(1e-6*ncsStats.sampleTSDiff);
    fid = fopen(fullfile(completeOutdir,'params.txt'),'w');
    fprintf(fid,'samplerate %d\n',samplerate);
    fclose(fid);
  end
  spikeClasses=cluster_class(:,1);
  spikeTimes=(cluster_class(:,2)*1000)+ncsStats.firstSampleTime; 
  
  clusters=setdiff(unique(spikeClasses),0);
  disp(clusters');
  if ~isempty(clusters)
    for cNum=1:length(clusters)
      % we have found another cell
      currentCell = currentCell+1;
      % need to add to events since we've found one cell
      if currentCell==1
	event_ms = getStructField(events,'mstime');
	pulse_micros = (b(1) +event_ms*b(2))*1000;
	pulse_sample = (pulse_micros-ncsStats.firstSampleTime)/ncsStats.sampleTSDiff;
	for e = 1:length(events)
	  events(e).unitoffset = pulse_sample(e);
	end
	events = replicateField(events,'unitfile',fullfile(completeOutdir,sprintf('unit_%s',outFileStem)));
      end
      % make an array with zeros for the whole length of recording    
      spike_present = zeros(1,dataLength);
      % add a 1 when that sample has a spike    
      cl=clusters(cNum);
      thisCellSpikes=spikeTimes(spikeClasses==cl);
      % subtract the starting time
      thisCellSpikes = (thisCellSpikes-ncsStats.firstSampleTime);
      % convert the spike times into samples 
      thisCellSamples = fix(thisCellSpikes/ncsStats.sampleTSDiff);
      spike_present(thisCellSamples) = 1;
      % save the data to file
      outFile = fullfile(completeOutdir,sprintf('unit_%s.%03d',outFileStem,c));
      fid = fopen(outFile,'wb','l');
      fwrite(fid,spike_present,'char');
      fclose(fid);
      % zip the file
      eval(['!gzip -f ' outFile]);
    end % loop over clusters for this channel
  end %the are clusters to loop over

  
  
end % loop over channels


