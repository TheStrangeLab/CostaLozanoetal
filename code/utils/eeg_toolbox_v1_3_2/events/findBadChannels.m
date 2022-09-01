function [bad,ranges,artBins] = findBadChannels(fileRoot,channels,binSizeMS,maxRange,minChanBad,maxBinBad,artifactChannels,artifactThresh)
%FINDBADCHANNELS - Identify bad channels of EEG data
%
% This function will detect bad channels from a session and both
% return and write those bad channels to fileRoot.bad_chan.  You
% can then use the result in a call to reref. 
%
% The basic method is to split each channel into bins and to check
% if the range of each bin exceeds maxRange, which would indicate a
% large fluctuation most likely due to a bad channel.  If you
% exceed the maxBinBad threshold for a percentage of bins being bad
% for a single channel, then that channel is marked bad.  If you
% provide artifactChannels and thresholds, then the function will
% ignore periods of eyeblinks in determining the bad channels.  If
% a bin is bad in more than minChanBad percent of the channels,
% then that bin will also be ignored.
% 
% FUNCTION:
%   bad = findBadChannels(fileRoot,channels,binSizeMS,maxRange,minChanBad,maxBinBad,artifactChannels,artifactThresh)
%
%
% INPUT ARGS:
%   fileRoot
%   channels
%   binSizeMS
%   maxRange
%   minChanBad
%   maxBinBad
%   artifactChannels
%   artifactThresh
%
% OUTPUT ARGS:
%   bad
%   ranges
%   artBins
%   
%
%

if ~exist('binSizeMS','var')
  binSizeMS = 1000;
end
if ~exist('maxRange','var')
  maxRange = 200;
end
if ~exist('minChanBad','var')
  minChanBad = .8;
end
if ~exist('maxBinBad','var')
  maxBinBad = .2;
end
if ~exist('artifactChannels','var')
  artifactChannels = [];
end
if ~exist('artifactThresh','var')
  artifactThresh = 100;
end


% make fake event
ev = struct('eegfile',{fileRoot},'eegoffset',{1});

% get data info
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(ev);

% get binsize
binsize = fix(binSizeMS*samplerate/1000);

% loop over channels files
fprintf('Channel: ')
for c = channels
  fprintf('%d',c);
  
  % load the data
  eeg = gete(c,ev,0,0,0);
  eeg = eeg{1};
  fprintf('.');
  
  % allocate matrix if first one
  if c == channels(1)
    % calc number of bins
    numbins = fix(length(eeg)/binsize);
    
    % allocate channel by bin matrix
    ranges = zeros(round(max(channels)),numbins);
  end
  
  % loop over bins
  % reshape for speed
  ranges(c,:) = range(reshape(eeg(1:numbins*binsize),binsize,numbins));
  %ranges(c,:) = kurtosis(reshape(eeg(1:numbins*binsize),binsize,numbins),0);
 
end
fprintf('Done\n');


% see if we are to get artifact info
if ~isempty(artifactChannels)
  blinks = addArtifacts(fileRoot,artifactChannels,artifactThresh);
  artInfo = zeros(1,length(eeg));
  artInfo(blinks) = 1;
  artBins = sum(reshape(artInfo(1:numbins*binsize),binsize,numbins))>0;
else
  artBins = zeros(1,numbins);
end

% set to 1/0
ctest = ranges;
ctest(ctest<maxRange) = 0;
ctest(ctest>0) = 1;

% do bad channels per bin percent checks
% this will ignore bins that are bad in almost all channels
ctest = ctest(:,find(mean(ctest,1)<minChanBad | ~artBins));
%ign = find(mean(ctest,1)>minChanBad);
%ctest(:,ign) = ones(size(ctest,1),length(ign))*0;

% do bad bins per channel percent checks to decide on bad channels
bad = find(mean(ctest,2)>maxBinBad);

% print out bad channels
fprintf('Bad Channels(%d): ',length(bad));
fprintf('%d ',bad);
fprintf('\n');

% write out bad channels
fid = fopen([fileRoot '.bad_chan'],'w');
fprintf(fid,'%d\n',bad);
fclose(fid);
