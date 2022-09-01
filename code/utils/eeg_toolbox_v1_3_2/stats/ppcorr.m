function ppcorr(dataroot,Channels,events,DurationMS,OffsetMS,BufferMS,varargin)
%PPCORR - Phase versus Power Correlation data preperation
%
% Calculate the data needed to calculate the correlation between 
% pre and post event power and
% phase reset.  This analysis helps to demonstrate that a phase reset
% that occurs due to an event is not due to evoked oscillations
% that increase power.
% 
% Power is measured as the t value of the paired ttest between pre
% and post event mean wavelet power, while phase is the mean
% Rayleigh statistic post event.
%
% FUNCTION:
%   ppcorr()
%
% INPUT VARS:
%
%
% OUTPUT VARS:
%
%
%


% make sure directory exists
[datadir,datafile] = fileparts(dataroot);
if ~exist(datadir,'dir')
  mkdir(datadir);
end
if ~exist('temp','dir')
  mkdir('temp');
end

% see if already done
if exist([dataroot '.mat'],'file')
  % already done
  disp([dataroot ' exists, so exiting...']);
  return
end

% default no merge
domerge = 0;
numskipped = 0;

% loop over each channel
for c = Channels
  % make sure not in progress
  filename = ['temp/' datafile '_' num2str(c) '.mat'];
  if exist(filename,'file')
    numskipped = numskipped + 1;
    continue;
  end
  
  % touch the file
  fid = fopen(filename,'wt');
  fclose(fid);

  % get power and phase for each event
  [phase,pow] = getphasepow(c,events,DurationMS,OffsetMS,BufferMS,varargin{:});
  
  % allocate space
  Ws = zeros(1,size(pow,2));
  Rbar = zeros(1,size(phase,2));
  
  % split in half
  midpoint = round(size(phase,3)/2);
  
  % calc mean power for both halves
  pre_pow = squeeze(mean(pow(:,:,1:midpoint),3));
  post_pow = squeeze(mean(pow(:,:,midpoint+1:end),3));
  
  % perform signtest, saving the sign
  for f = 1:size(pow,2)
    % do sign test
    [p_temp,h_temp,stats] = signrank_ci(post_pow(:,f),pre_pow(:,f));
    
    % save the sign and direction
    Ws(f) = stats.wdiff;
  end
  
  % use second half of phase and calc the rayleigh stat
  for f = 1:size(phase,2)
    % get the rayleigh for each point in time
    [p,rb_temp] = rayleigh(squeeze(phase(:,f,midpoint+1:end)));
    
    % save mean of Rbar
    Rbar(f) = mean(rb_temp);
  end

  % save the channel to a file
  save(filename,'Ws','Rbar');
  
  % see if merge
  if c == Channels(end)
    domerge = 1;
  end
end % channels

% see if we merge the data
if (domerge == 1 | numskipped == length(Channels))
  % doing merge
  fprintf('\nMerging data...\n');

  % loop over the channels
  fprintf('Real data (%d): ',max(Channels));
  for c = Channels
    fprintf('%d ',c);
    % load the file
    filename = ['temp/' datafile '_' num2str(c) '.mat'];
    load(filename)

    % allocate space
    if ~exist('Ws_all','var')
      Ws_all = zeros(max(Channels),length(Ws));
      Rbar_all = zeros(max(Channels),length(Rbar));
    end
    
    % combine the data
    Ws_all(c,:) = Ws;
    Rbar_all(c,:) = Rbar;
  end

  % save the combined data
  save(dataroot,'Ws_all','Rbar_all');
  
  % remove the temp files
  delete(['temp/' datafile '*']);

end






