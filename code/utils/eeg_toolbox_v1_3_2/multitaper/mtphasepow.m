function [varargout] = mtphasepow(chan,events,DurationMS,OffsetMS,BufferMS,varargin)
%MTPHASEPOW - Calculate wavelet phase and power for a set of events.
%
% Calculate multitaper phase and power as a function of time and frequency for a
% single electrode.
%
% FUNCTION: 
%   [phase,pow,kInd] = mtphasepow(chan,events,DurationMS,OffsetMS,BufferMS,varargin)
%
% INPUT ARGs:
%   chan = 2;
%   events = events;
%   DurationMS = 2000;
%   OffsetMS = 0;
%   BufferMS = 1000;
%
%   OPTIONAL PARAMS:
%     'freqs'
%     'bandwidth'
%     'windowsize'
%     'filtfreq'
%     'filttype'
%     'filtorder'
%     'resampledrate' - resample applied before calculating phase/power
%     'downsample' - decimate applied after calculating power
%     'powonly'
%     'usesingles'
%     'kthresh'- Kurtosis threshold to throw out events.
%
% OUTPUT ARGS:
%   phase- (Events,Freqs,Time)
%   pow- (Events,Freqs,Time)
%   kInd - Logical indexes of the events not thrown out due to kurtosis
%

%
% Changes:
%
% 9/15/05 - PBS - Added downsampling via decimate following power calculation.
% 9/15/05 - PBS - Return the logical index of the events not thrown out
%                 with the kurtosis thresh.


% set the defaults
freqs = eeganalparams('freqs');
bandwidth = eeganalparams('bandwidth');
windowsize = eeganalparams('windowsize'); 
filtfreq =  [];
filttype = 'stop';
filtorder = 4;
resampledrate = [];
dsample = [];
powonly = 0;
usesingles = 0;
kthresh = [];

% process the varargs
i = 1;
while length(varargin)> 0 & i<=length(varargin)
  switch lower(varargin{i})
   case 'freqs'
    i = i+1;
    freqs = varargin{i};
    i = i+1;
   case 'bandwidth'
    i = i+1;
    bandwidth = varargin{i};
    i = i+1;
   case 'windowsize'
    i = i+1;
    windowsize = varargin{i};
    i = i+1;
   case 'filtfreq'
    i = i+1;
    filtfreq = varargin{i};
    i = i+1;
   case 'filttype'
    i = i+1;
    filttype = varargin{i};
    i = i+1;
   case 'filtorder'
    i = i+1;
    filtorder = varargin{i};
    i = i+1;
   case 'resampledrate'
    i = i+1;
    resampledrate = varargin{i};
    i = i+1;
   case 'downsample'
    i = i+1;
    dsample = varargin{i};
    i = i+1;
   case 'powonly'
    i = i+1;
    powonly = 1;
   case 'usesingles'
    i = i+1;
    usesingles = 1;
   case 'kthresh'
    i = i+1;
    kthresh = varargin{i};
    i = i+1;
    
   otherwise
    error(['Error processing vararg: ' num2str(i)]);
  end
end

% get some parameters
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(1));
rate = eegparams('samplerate',fileparts(events(1).eegfile));
samplerate = rate;
if ~isempty(resampledrate)
  resampledrate = round(resampledrate);
  samplerate = resampledrate;
  rate = resampledrate;
end

% convert the durations to samples
%duration = fix((DurationMS+(2*BufferMS))*rate/1000);
%offset = fix((OffsetMS-BufferMS)*rate/1000);
buffer = fix((BufferMS)*rate/1000);

% load the eeg
eeg = gete_ms(chan,events,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,filtfreq,filttype,filtorder,samplerate);

% see if throw out events with weird kurtosis
if ~isempty(kthresh)
  startsize = size(eeg,1);
  k = kurtosis(eeg(:,buffer+1:end-buffer)');
  goodInd = k<=kthresh;
  kInd = goodInd;
  eeg = eeg(goodInd,:);
  sizediff = startsize - size(eeg,1);
  if sizediff > 0
    fprintf('Threw out %d events due to kurtosis...\n',sizediff);
  end
else
  kInd = logical(ones(size(eeg,1),1));
end

% zero out the data
if usesingles
  if ~powonly
    phase = single(zeros(size(eeg,1),size(freqs,2),size(eeg,2)));
  end
  pow = single(zeros(size(eeg,1),size(freqs,2),size(eeg,2)));
else
  if ~powonly
    phase = zeros(size(eeg,1),size(freqs,2),size(eeg,2));
  end
  pow = zeros(size(eeg,1),size(freqs,2),size(eeg,2));
end
  
% get the power
if length(events) > 1
  fprintf('Events: %g\n',size(eeg,1)); 
end
work = [];
for j = 1:size(eeg,1)
  if powonly
    [pow(j,:,:),phase,work] = mtenergyvec(eeg(j,:),freqs,samplerate,bandwidth,windowsize,work);
  else
    [pow(j,:,:),phase(j,:,:),work] = mtenergyvec(eeg(j,:),freqs,samplerate,bandwidth,windowsize,work);
  end
  
  if length(events) > 1
    fprintf('%g ',j);
  end
end

if length(events) > 1
  fprintf('\n');
end

% remove the buffer
pow = pow(:,:,buffer+1:end-buffer);

% see if decimate power
if ~isempty(dsample)
  % set the downsampled duration
  dmate = round(samplerate/dsample);
  dsDur = ceil(size(pow,3)/dmate);

  if usesingles
    precision = 'single';
  else
    precision = 'double';
  end
  dpow = zeros(size(pow,1),size(pow,2),dsDur,precision);

  % Must log transform before modifying
  pow(pow<=0) = eps;
  pow = log10(pow);
  
  % loop and decimate
  for e = 1:size(pow,1)
    for f = 1:size(pow,2)
      dpow(e,f,:) = decimate(double(pow(e,f,:)),dmate);
    end
  end
  
  % replace old pow with new
  pow = dpow;
  clear dpow;
  
  % convert back to no-log
  pow = 10.^pow;
  
end


if ~powonly
  phase = phase(:,:,buffer+1:end-buffer);
  varargout(1) = {phase};
  varargout(2) = {pow};
  varargout(3) = {kInd};
else
  varargout(1) = {pow};
  varargout(2) = {kInd};
end

