function [m,s] = calcNormFact(chan,events,num,DurationMS,filtFreq)
%CALCNORMFACT - Calculate a mean and std for power Z-transform
%
%
% FUNCTION:
%   [m,s] = calcNormFact(chan,events,num,DurationMS,filtFreq)
%
%
%
%

if nargin < 5
  filtFreq = [];
end

% get samplerate info
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(1));
duration = round(DurationMS*samplerate/1000);

freqs = eeganalparams('freqs');

% get unique file names
uFiles = unique(getStructField(events,'eegfile',''));

% get start and end for each filename
epochs = [];
for f = 1:length(uFiles)
  % get all offsets for each unique file
  offsets = getStructField(events,'eegoffset',['strcmp(eegfile,''' uFiles{f} ''')']);
  
  % append epochs from file
  fepochs = min(offsets):duration:max(offsets);
  epochs = [epochs [ones(1,length(fepochs)).*f;fepochs]];
end

% randomize which ones to use
poss_epochs = size(epochs,2);
fprintf('Possible Epochs: %d\n',poss_epochs);

if isempty(num) | num == 0
  num = poss_epochs
end

x = randperm(poss_epochs);
epochs = epochs(:,x(1:num));

% create the events
epoch_events = struct('eegfile',{uFiles{epochs(1,:)}},'eegoffset',mat2cell(epochs(2,:)));
  
% loop over epochs, calculating power
m = zeros(length(freqs),1);
s = zeros(length(freqs),1);
fprintf('Epochs(%d): ',num);
for e = 1:length(epoch_events)
  fprintf('%d ',e);
  % get the power
  [phase,pow] = getphasepow(chan,epoch_events(e),DurationMS,0,1000,'filtfreq',filtFreq);
  
  mpow = mean(squeeze(pow),2);
  spow = std(squeeze(pow),0,2);
  m = m + mpow;
  s = s + spow;
end

fprintf('\n');

% get average from sums
m = m/num;
s = s/num;




% calc mean and std across events and time
%m = mean(mean(pow,3),1);
%s = std(std(pow,0,3),0,1);
