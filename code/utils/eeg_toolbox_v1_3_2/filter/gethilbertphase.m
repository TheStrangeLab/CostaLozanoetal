function [phase,amp,eeg] = gethilbertphase(chan,events,DurationMS,OffsetMS,BufferMS,bandpassFreq,notchFilterFreq)
%GETHILBERTPHASE - uses a hilbert transformation to get the instantaneous phase
%in a particular bandpassed frequency band.
%
% FUNCTION: 
%function phase = gethilbertphase(chan,events,DurationMS,OffsetMS,BufferMS,bandpassFreq,notchFilterFreq)
%
% INPUT ARGs:
%   chan = 2;
%   events = events;
%   DurationMS = 2000;
%   OffsetMS = 0;
%   BufferMS = 1000;
%   bandpassFreq=[4 8]
%   notchFilterFreq=60
%
% OUTPUT ARGS:
%   phase- (Events,Time)


[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(1));
rate = eegparams('samplerate',fileparts(events(1).eegfile));
%rate=500;

% convert the durations to samples
duration = fix((DurationMS+(2*BufferMS))*rate/1000);
offset = fix((OffsetMS-BufferMS)*rate/1000);
buffer = fix((BufferMS)*rate/1000);
 

% load the eeg
eeg = gete_ms(chan,events,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0, ...
                bandpassFreq, 'bandpass',2,rate);
if iscell(eeg)
  if length(eeg)~=1
    error(['If you are asking for the entire file''s duration, you must only ' ...
           'pass in one event, otherwise this code doesn''t support it.']);
  end
  eeg=eeg{1};
end

for row=1:size(eeg,1)
  if exist('notchFilterFreq')
%    eeg(row,:)=sineFitNotch(eeg(row,:),samplerate,notchFilterFreq);
    eeg(row,:)=buttfilt(eeg(row,:),[-1 1]+notchFilterFreq,samplerate,'stop',2);
  end
  h=hilbert(eeg(row,:));
  phase(row,:)=angle(h);
  amp(row,:)=abs(h);
end

%remove the buffering
phase = phase(:,buffer+1:end-buffer);
amp = amp(:,buffer+1:end-buffer);

