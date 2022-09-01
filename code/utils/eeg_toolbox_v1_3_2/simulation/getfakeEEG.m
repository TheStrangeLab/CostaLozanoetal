function eeg=getfakeEEG(numEvents,durationMS,maxAmp,samplingRate)
%GETFAKEEEG Create simulated EEG event data.
%
% This function creates fake, random background EEG to be used for
% simulations; Every  epoch has a different random phase and
% frequency. The amplitude of the EEG makes sure its power follows
% a 1/f distribution.
%
% FUNCTION:
%   eeg=getfakeEEG(numEvents,durationMS,maxAmp,samplingRate)
%
% INPUT ARGS:
%   numEvents = 10    % number of epochs to create
%   duration = 2000   % duration of the event in ms
%   maxAmp = 0.5      % this determines the scale of the EEG. To
%      approximate scalp EEG a value of 0.5 works reasonably well,
%      for intracranial most likely a slightly larger value is
%      required. The number is the amplitude of a composite sine wave
%      at 0.1 Hz (which will fall off with frequency to create a 1/f
%      spectrum)
%   samplingRate=256  % the samplingrate at which the signal
%      should be collected
%
% OUTPUT ARGS:
%    eeg (numEvents,samples)
%
%

%initialize the random number generator
rand('state',sum(100*clock));

if nargin < 4
  samplingRate = 256;
  if nargin < 3
    maxAmp = 0.5
  end
end

% set some variables to control the gake eeg
maxFreq = 100; %Hz, the maximum frequency of the background
               %spectrum
minFreq = 0.1;
numSineWaves = 500; %number of since waves used in computing the
                   %background spectrum

% create the time vector
duration = fix(durationMS*samplingRate/1000);
eeg = zeros(numEvents,duration);
timeVec = linspace(0,durationMS/1000,duration);
timeVec = timeVec(ones(numEvents,1),:);

% Set the amp factor.
% The factor with which each 1/f will be multiplied to yield the
% correct amplitude (i.e., the slope of the 1/f line)
ampfactor = (maxAmp/minFreq)*ones(numEvents,1); 

% create the background by summing 50 sinusoids of random frequency
% and phase

% pick the random frequencies and phase for all waves going into
% each event
freq = abs((maxFreq-minFreq)*rand(numEvents,numSineWaves) + minFreq);
rndPhase = 2*pi*rand(numEvents,numSineWaves);

% allocate space for the eeg
eeg = zeros(numEvents,duration);

% loop over waves to add
for n = 1:numSineWaves
  
  % determine the amp factor
  tFactor = ampfactor./sqrt(freq(:,n));
  tempFactor = tFactor(:,ones(duration,1));
  
  % determine the sine waves
  tWave = 2*pi*freq(:,n);
  tPhase = rndPhase(:,n);
  tempWave = sin(tWave(:,ones(duration,1)).*timeVec + tPhase(:,ones(duration,1)));
  
  % add it to the eeg
  eeg = eeg + tempFactor.*tempWave;
  
end
