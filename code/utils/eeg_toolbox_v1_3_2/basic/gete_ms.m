function EEG=gete_ms(channel,events,DurationMS,OffsetMS,BufferMS,filtfreq,filttype,filtorder,resampledRate, RelativeMS)
%GETE_MS - Get EEG event data based on MSec ranges instead of samples.
% 
% Returns data from an eeg file.  User specifies the channel,
% duration, and offset along with an event.  The event struct MUST
% contain both 'eegfile' and 'eegoffset' members.
%
% You can optionally resample the data with the Signal Toolbox's
% resample function.  The resampling occurs following the filter.
% Use the resampling with caution because the we have not provent
% that the spectral properties of the data do not change dramatically.
%
% NOTE: All events must have the same sampling rate.
%
% FUNCTION:
%   EEG=gete_ms(channel,events,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate,RelativeMS)
%
% INPUT ARGS:
%   channel = 3;            % the electrode #
%   events = events(8:12);  % event struct to extract [eegfile eegoffset]
%   durationMS = 2000;      % signal time length in samples
%   offsetMS = 0;           % offset at which to start in samples
%   bufferMS = 1000;        % buffer (needed for filtering or resampling)
%                           %   default is 0
%   filtfreq = [58 62];     % Filter freq (depends on type, see buttfilt)
%                           %   default is []
%   filttype = 'stop';      % Filter type (see buttfilt)
%   filtorder = 1;          % Filter order (see buttfilt)
%   resampledRate = 200;    % Sample rate of the returned data
%   RelativeMS = [-200 0];  % Range for use with the relative subtraction
%
% OUTPUT ARGS:
%   EEG(Trials,Time) - The data from the file
%

% 12/18/07 - MvV - changed the indices into readbytes when saving
% to EEG, such that it was always fit and not be affected by
% rounding differences.
% 11/29/04 - PBS - Changed round to fix to fix range problem
% 4/20/04 - PBS - Added Relative Range subtraction

% check the arg
if nargin < 10
  RelativeMS = [];
  if nargin < 9
    resampledRate = [];
    if nargin < 8
      filtorder = 1;
      if nargin < 7
	filttype = 'stop';
	if nargin < 6
	  filtfreq = [];
	  if nargin < 5
	    buffer = 0;
	    if nargin<4 
	      offset=0; 
	    end; 
	  end
	end
      end
    end
  end
end

% get initial data info
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(1));
samplerate = round(samplerate);
if isempty(resampledRate)
  resampledRate = samplerate;
end
resampledRate = round(resampledRate);

% base final datasize on resampled data
final_duration = fix((DurationMS)*resampledRate/1000);
final_offset = fix((OffsetMS)*resampledRate/1000);
final_buffer = fix((BufferMS)*resampledRate/1000);

% allocate space
EEG = zeros(length(events),final_duration);

for e = 1:length(events)
  
  % get data info
  [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(e));
  samplerate = round(samplerate);
  
  % set event durations from rate
  duration = fix((DurationMS+(2*BufferMS))*samplerate/1000);
  offset = fix((OffsetMS-BufferMS)*samplerate/1000);
  buffer = fix((BufferMS)*samplerate/1000);

  % set the channel filename
  eegfname=sprintf('%s.%03i',events(e).eegfile,channel);

  eegfile=fopen(eegfname,'r','l'); % NOTE: the 'l' means that it
                                   % came from a PC! 
  if(eegfile== -1)
    eegfname=sprintf('%s%03i',events(e).eegfile,channel); % now try unpadded lead#
    eegfile=fopen(eegfname,'r','l');
  end
  if(eegfile==-1)
    eegfname=sprintf('%s.%i',events(e).eegfile,channel); % now try unpadded lead#
    eegfile=fopen(eegfname,'r','l');
  end

  % tell if not open
  if eegfile==-1
    % did not open
    error('ERROR: EEG File not found for event(%d): %s.\n',e,events(e).eegfile);
  end
    
  % read the eeg data
  thetime=offset+events(e).eegoffset;
  fseek(eegfile,nBytes*thetime,-1);

  readbytes=fread(eegfile,duration,dataformat)';
  if length(readbytes)~=duration
    warning([eegfname ' only ' num2str(length(readbytes)) ' of ' num2str(duration) ' samples read for event ' num2str(e) ', appending zeros']);
    readbytes=[readbytes zeros(1,duration-length(readbytes))];
  end
  
  
  % close the file
  fclose(eegfile);

  % see if filter the data
  if length(filtfreq) > 0
    readbytes=buttfilt(readbytes,filtfreq,samplerate,filttype,filtorder);
  end

  % see if resample
  if resampledRate ~= samplerate
    % do the resample
    readbytes = resample(readbytes,round(resampledRate),round(samplerate));
    %readbytes = dsamp(readbytes,round(samplerate),round(resampledRate));
  end
  
  % append the data
  EEG(e,:) = readbytes(final_buffer+1:(final_buffer+final_duration));
end

% see if take relative baseline correction
if length(RelativeMS) == 2
  % get the average for the range
  relative = fix((RelativeMS)*resampledRate/1000);
  relative = relative - final_offset + 1;
  relative(2) = relative(2) - 1;
  
  % calculate the relative
  releeg = mean(EEG(:,(relative(1):relative(2))),2);

  % subtract the baseline
  for e = 1:size(EEG,1)
    EEG(e,:) = EEG(e,:)-releeg(e);     
  end
end 

% do the gain multiplication
EEG=EEG.*gain;

