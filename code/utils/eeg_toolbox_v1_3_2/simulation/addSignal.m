function neweeg =  addSignal(eeg,frequencies,amp,signalTime,sampleRate,startPhase)
% ADDSIGNAL - add a signal to the fake background EEG
%
% this function adds a desired signal to a (fake) piece of
% EEG,e.g. a piece of background EEG that has been generated using
% getFakeEEG
%
% FUNCTION:
% function neweeg =  addSignal(eeg,frequencies,amp,signalTime,sampleRate,startPhase)
%
% INPUT ARGs
%     eeg = a stretch of EEG signal, as resulting from either
%     getFakeEEG or gete, which has as its rows the different
%     events, and as its columns time
%     frequencies  = [10]          - frequencies at which to add a signal
%     amp = [1]                - muV; amplitude to give the signal
%     signalTime = [100 400] - ms the time in ms during the trial during
%     which the signal takes place
%     sampleRate = 256       - the samplerate of the signal
%     startPhase = [0]         - a number between 0 and 2*pi that
%     indicates a which phase the signal should start (optional)
% OUTPUT ARGs
%     neweeg = the input signal with the phasic EEG peak as
%     specified by the input arguments added to it
 
if nargin < 6
   startPhase = zeros(1,length(frequencies));
 end
 
numEvents = size(eeg,1);
epochLength = size(eeg,2)*(1/sampleRate);
timeVec = (1/sampleRate):(1/sampleRate):epochLength;
signal = zeros(size(eeg));
neweeg = eeg; %copy the background signal to the new variable

for n = 1:length(frequencies)
  timeInd = fix((signalTime(n,:).*sampleRate)./1000);
  timeInd(timeInd==0) = 1; % make sure the index does not go to 0
  baseSig = amp(n)*sin((2*pi*frequencies(n)).*timeVec(timeInd(1):timeInd(2)) + startPhase(n));
  signal(:,timeInd(1):timeInd(2)) = baseSig(ones(numEvents,1),:);
  neweeg = neweeg + signal;
end


