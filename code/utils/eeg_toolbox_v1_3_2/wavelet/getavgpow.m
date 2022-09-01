function pow = getavgpow(chan,events,DurationMS,OffsetMS,BufferMS,varargin)
%getavgpow - Calculate time-averaged power for a window around each
%event.
%
% FUNCTION: 
%   [pow,allpows] = getavgpow(chan,events,DurationMS,OffsetMS,BufferMS,varargin)
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
%     'width'
%     'filtfreq'
%     'filttype'
%     'filtorder'
%     'averageAcrossFreqs'
%
% OUTPUT ARGS:
%   pow- (Events,Freqs,Time)
%   allpows-



% set the defaults
freqs = eeganalparams('freqs');
width = eeganalparams('width');
filtfreq =  [];
filttype = 'stop';
filtorder = 1;
averageAcrossFreqs=0;        

% process the varargs
i = 1;
while length(varargin)> 0 & i<=length(varargin)
  switch lower(varargin{i})
   case 'freqs'
    i = i+1;
    freqs = varargin{i};
    i = i+1;
   case 'width'
    i = i+1;
    width = varargin{i};
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
   case 'averageacrossfreqs'
    i = i+1;
    averageAcrossFreqs=1;        
   otherwise
    error(['Error processing vararg: ' num2str(i)]);
  end
end

% get some parameters
rate = eegparams('samplerate',fileparts(events(1).eegfile));
samplerate = rate;

% convert the durations to samples
duration = round((DurationMS+(2*BufferMS))*rate/1000);
offset = round((OffsetMS-BufferMS)*rate/1000);
buffer = round((BufferMS)*rate/1000);



% load the eeg
eeg = gete(chan,events,duration,offset,buffer,filtfreq,filttype,filtorder);

% zero out the data
pow = zeros(size(eeg,1),size(freqs,2));

% get the power
fprintf('Events: %g\n',size(eeg,1)); 
for j = 1:size(eeg,1)
  [phase,p] = multiphasevec2(freqs,eeg(j,:),samplerate, ...
                                     width);
  %return the average of the non buffered time

  pow(j,:)=mean(p(:,buffer+1:end-buffer),2)';
  fprintf('%g ',j);
end

fprintf('\n');

if averageAcrossFreqs
  pow=mean(pow,2);
end


