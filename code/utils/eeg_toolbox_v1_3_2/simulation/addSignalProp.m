function neweeg =  addSignalProp(eeg,frequency,ampFr,signalTime,sampleRate,startPhase)
% this function adds a desired signal to a (fake) piece of
% EEG,e.g. a piece of background EEG that has been generated using
% getFakeEEG. The amplitude is taken as a fraction of the generated EEG
%
% INPUT ARGs
%     eeg = a stretch of EEG signal, as resulting from either
%     getFakeEEG or gete, which has as its rows the different
%     events, and as its columns time
%     frequency  = 10          - frequency at which to add a signal
%     ampFr = [0.01]                - amplitude to give the signal as
%     a fraction of the background at that frequency; this can be
%     either a scalar, in which the same amplitude will be applied
%     to every piece of EEG (every event) or a vector which
%     contains a different amplitude for every event
%     signalTime = [100 400] - ms the time in ms during the trial during
%     which the signal takes place; this can also have an extra (second)
%     dimension which tells us how long the signal should be for
%     every event
%     sampleRate = 256       - the samplerate of the signal
%     startPhase = 0        - a number between 0 and 2*pi that
%     indicates a which phase the signal should start (optional)
% OUTPUT ARGs
%     neweeg = the input signal with the phasic EEG peak as
%     specified by the input arguments added to it
 
if nargin < 6
   startPhase = 0;
 end
 
numEvents = size(eeg,1);
epochLength = size(eeg,2)*(1/sampleRate);
timeVec = (1/sampleRate):(1/sampleRate):epochLength;
signal = zeros(size(eeg));
neweeg = eeg; %copy the background signal to the new variable

  % determine the actual amplitude for this frequency (using the
  % way the background EEG was constructed, i.e., 50 sinusoids that
  % have the amplitude fall-off described by the function below)
  amplBackgr = ((.5/.1)*50./sqrt(frequency));
  if size(signalTime,2)==1 %we don't have different timecourse for
                        %different events
    timeInd = fix((signalTime.*sampleRate)./1000);
  else
    timeInd = squeeze(fix((signalTime.*sampleRate)./1000));
  end
  timeInd(timeInd==0) = 1; % make sure the index does not go to 0
  if size(ampFr,1)==1 % one amplitude specified for all events
    amp = ampFr*amplBackgr;
    if size(signalTime,1)==1
      signal(:,timeInd(1):timeInd(2)) = repmat(amp*sin((2*pi*frequency+ startPhase).*timeVec(timeInd(1):timeInd(2))),numEvents,1);
    else
      for e = 1:numEvents
	signal(e,timeInd(e,1):timeInd(e,2)) = amp*sin((2*pi*frequency+ startPhase).*timeVec(timeInd(e,1):timeInd(e,2)));
      end
    end
  else %different amplitudes for different events
    amp = ampFr*amplBackgr;
    if size(signalTime,1)==1
      signal(:,timeInd(1):timeInd(2)) = repmat(amp,1,length(timeInd(1):timeInd(2))).*repmat(sin((2*pi*frequency+ startPhase).*timeVec(timeInd(1):timeInd(2))),numEvents,1);    
    else
      for e = 1:numEvents
           signal(e,timeInd(e,1):timeInd(e,2)) = amp(e).*sin((2*pi*frequency+ startPhase).*timeVec(timeInd(e,1):timeInd(e,2))); 
      end
    end
  end
  neweeg = neweeg + signal;




