function [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event)
%GETRATEANDFORMAT - Get the samplerate, gain, and format of eeg data.
%
% function [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event)
%

if isstr(event)
  paramdir = event;
else
  paramdir = fileparts(event.eegfile);
end

samplerate = eegparams('samplerate',paramdir);
gain = eegparams('gain',paramdir);

dataformat='short';
tmpdataformat= eegparams('dataformat',paramdir);
if ~isempty(tmpdataformat)
  dataformat=tmpdataformat;
end

switch dataformat
 case 'single'
  nBytes=4;
 case {'short','int16'}
  nBytes=2;
 case 'double'
  nBytes=8;
 otherwise
  error('BAD DATA FORMAT!!!!!');
end
