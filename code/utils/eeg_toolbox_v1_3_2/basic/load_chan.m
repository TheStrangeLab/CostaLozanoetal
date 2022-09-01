function [eeg,stats] = load_chan(fileroot,channel,OffsetMS,DurationMS,varargin)
	%
	%LOAD_CHAN   Load eeg data from one channel.
	%   EEG = LOAD_CHAN(FILEROOT,CHANNEL,OFFSETMS,DURATIONMS,VARARGIN)
	%   loads the file
	
	% process the input arguments
	pnames = {'relativems', 'resampledrate', 'bufferms', 'filttype', 'filtfreq', 'filtorder'};
	dflts = {[], 500, 1000, '', [], []};
	[eid,emsg,RelativeMS,resampledRate,BufferMS,filttype,filtfreq,filtorder] = getargs(pnames, dflts, varargin{:});
	
	[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(fileparts(fileroot));
	
	% set event durations from rate
  duration = fix(DurationMS*samplerate/1000);
  offset = fix((OffsetMS-BufferMS)*samplerate/1000);
  buffer = fix((BufferMS)*samplerate/1000);
  totduration = duration + buffer*2;
	relative = fix(OffsetMS+RelativeMS.*samplerate/1000);

  % set the channel filename
  eegfname = sprintf('%s.%03i',fileroot,channel);
  fid = fopen(eegfname,'r','l');

  if fid==-1
	  % now try no dot before lead #
    eegfname = sprintf('%s%03i',fileroot,channel); 
    fid = fopen(eegfname,'r','l');
  end
  if fid==-1
	  % now try unpadded lead #
    eegfname = sprintf('%s.%i',fileroot,channel); 
    fid = fopen(eegfname,'r','l');
  end
  if fid==-1
    % did not open
    error('ERROR: EEG File not found: %s.\n', fileroot);
  end
  
  % find the starting place
  fseek(fid,nBytes*offset,-1);

	% read the data
  readbytes = fread(fid,totduration,dataformat)';
  if length(readbytes)~=totduration
		warning('%s: only %d of %d samples read; appending zeros', ...
		eegfname, length(readbytes), totduration);

    readbytes = [readbytes zeros(1,totduration-length(readbytes))];
  end
  
  % close the file
  fclose(fid);

  % see if filter the data
  if ~isempty(filtfreq)
    readbytes = buttfilt(readbytes,filtfreq,samplerate,filttype,filtorder);
  end

  % see if resample
  if resampledRate ~= samplerate
    % do the resample
    readbytes = resample(readbytes,round(resampledRate),round(samplerate));
    %readbytes = dsamp(readbytes,round(samplerate),round(resampledRate));

		offset = fix((OffsetMS-BufferMS)*resampledRate/1000);
		duration = fix((DurationMS)*resampledRate/1000);
		buffer = fix((BufferMS)*resampledRate/1000);
		relative = fix(OffsetMS+RelativeMS*resampledRate/1000);
		samplerate = resampledRate;
  end
  
  % remove the buffer
  eeg = readbytes(buffer+1:(buffer+duration));
	offset = fix((OffsetMS-BufferMS)*samplerate/1000);
	% see if take relative baseline correction
	if length(RelativeMS) == 2
	  % get the average for the range
	  relative = relative - offset + 1;
	  relative(2) = relative(2) - 1;

	  % calculate the relative
	  releeg = mean(eeg(relative(1):relative(2)));

		% subtract the baseline
		eeg = eeg - releeg;
	end

	% do the gain multiplication
	eeg = eeg.*gain;
	
	% pass out a few things so we don't have to recalculate
	if nargout==2
		stats = struct('samplerate', samplerate, 'offsetMS', OffsetMS, 'durationMS', DurationMS);
	end
	