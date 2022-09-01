function newdata = removeBlinks(data,eogChannels,samplerate,winsizeMS,overlapMS,corrThresh)
%REMOVEBLINKS - Remove eyeblinks from continuous data using PCA
%
% Requires EEG Lab Toolbox for PCA functions.
%
% FUNCTION:
%   newdata = removeBlinks(data,eogChannels,samplerate,winsizeMS,overlapMS,corrThresh)
%
% INPUT ARGS:
%
%
% OUTPUT ARGS:
%
%


% set the default values
if ~exist('eogChannels','var') | isempty(eogChannels)
  eogChannels = [8,26,125,128];
end
if ~exist('samplerate','var') | isempty(samplerate)
  samplerate = 500;
end
if ~exist('winsizeMS','var') | isempty(winsizeMS)
  winsizeMS = 2000;
end
if ~exist('overlapMS','var') | isempty(overlapMS)
  overlapMS = 500;
end
if ~exist('corrThresh','var') | isempty(corrThresh)
  corrThresh = .4;
end

% convert MS to samples for winsize and overlap
winsize = fix(winsizeMS*samplerate/1000);
overlap = fix(overlapMS*samplerate/1000);

% allocate for new data
newdata = zeros(size(data),'single');

% set up tapers
taper = [0:1/(overlap-1):1 ones(1,(winsize-2*overlap)) 1:-1/(overlap-1):0];
taper = taper(ones(size(data,1),1),:);
start_taper = [ones(1,overlap) ones(1,(winsize-2*overlap)) 1:-1/(overlap-1):0];
start_taper = start_taper(ones(size(data,1),1),:);

% loop over window, combining overlaps
startInd = 1;
endInd = startInd + winsize-1;
while endInd <= size(data,2)
  % calc the pca of the segment
  [pc,eigenvec,sv]=runpca(data(:,startInd:endInd)');
  
  % calc the correlations with the EOG chans
  toRemove = [];
  rho = corr(data(eogChannels,startInd:endInd)',eigenvec).^2;
  [x,toRemove] = find(rho >= corrThresh);
  toRemove = unique(toRemove);
 
  % remove components if necessary
  if ~isempty(toRemove)
    % set those to zero
    eigenvec(:,toRemove) = 0;
    
    % project back
    tdata = [eigenvec*pc]';

    fprintf('-');
  else
    % just set back to original data
    tdata = data(:,startInd:endInd);

    fprintf('%g',max(max(rho)));
  end
    
  %keyboard
  
  % add in data to correct spot with a nice taper
  if startInd == 1
    % is first time, so don't taper on the left
    newdata(:,startInd:endInd) = newdata(:,startInd:endInd) + tdata.*start_taper;
  else
    % taper on both sides
    newdata(:,startInd:endInd) = newdata(:,startInd:endInd) + tdata.*taper;
  end
  
  % advance the window
  startInd = endInd - overlap + 1;
  lastEnd = endInd;
  endInd = startInd + winsize - 1;

  %[pc,eigenvec,sv]=runpca(squeeze(art_eeg(:,nc,:))');
  %eigenvec(:,correlatingComps) = 0;
  %ewdata = eigenvec*pc;
  
end

% see if there's any data left
if lastEnd < size(data,2)
  % append last data chunk unanalyzed
  % must taper it
  startInd = lastEnd - overlap + 1;
  endInd = size(data,2);
  duration = length(startInd:endInd);
  taper = [0:1/(overlap-1):1 ones(1,(duration-overlap))];
  taper = taper(ones(size(data,1),1),:);

  % Add in the data 
  tdata = data(:,startInd:endInd);
  newdata(:,startInd:endInd) = newdata(:,startInd:endInd) + tdata.*taper;
end

fprintf('\n');


