function calcPepisode_old(events,chan,outdir,filtFreq,freqs)
%CALCPEPISODE - Calculate Pepisode on entire EEG files.
%
% This function goes through all the EEG files that are specified by the input events structure and calculates the Pepisode transform of that. The output is a binary file which for each point in time contains a 0 if the power in that particular frequency is lower than baseline and a 1 if the power of that frequency is higher than the baseline. In order to get a P value for a certain time segment, use the getp_ms function
%
% The P value can be deduced from the 0 and 1 vector using
% u = find(vector);
% if (isempty(u)), P = 0; else P = length(u)/(length(vector); end
%
%
%   The output is the file root with 'union' added.
%    Channels is a vector with the desired channels to analize
%   The script uses by default the 95% level for the threshold. It
%   also uses the chi_squarefit function to fit the wavelet power
%   spectrum. If another type of fit is desired, this can be done
%   by simply replacing this function call with another function
%   call that produces a similar output.
%
%
%
% FUNCTION:
%   calcPepisode(events,chan,outdir,filtFreq,freqs)
%
% INPUT ARGS:
%   events = events - events that need a union vector
%   chan = [23 25]- the channel(s) to compute
%   outdir = 'pepisode_standard' - the directory in which to place
%   the output, if it doesn't contain an absolute path, it will be
%   placed in the designated directory relative to the current directory
%   filtFreq = [58 62]  = frequencies at which to apply a notch filter
%   freqs = (2^(1/8)).^(8:42) = frequencies at which union vectors
%   will be computed (SHOULD BE LOGARITHMICALLY SPACED)
%


if nargin < 5
  freqs = eeganalparams('freqs');
  if nargin < 4
    filtFreq = [];
  end
end

% set some vars
width = eeganalparams('width');

% get data info
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(1));

shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

outdir = fullfile('pepisode',outdir);

% make sure the required output directories exist, if they do not exist, create them
if ~exist(outdir,'dir')
  mkdir(outdir)
end

% get the filenames
allFile = getStructField(events,'eegfile','');
efnames = unique(allFile);
if isempty(efnames{1})
    efnames = efnames(2:end);
end

% see if already in progress
% get the file parts
[fpath,fname,fext] = fileparts(efnames{1});

% set the extension and filename
fext = sprintf('.%03i',chan);
filename = fullfile(outdir,[fname '_union' fext]);

if exist(filename,'file')
  % already done
  return
end

% load all the unique file data
fprintf('Loading EEG...');
t0 = clock;
eeg = gete(chan,events,0,0,0,filtFreq,'stop',1);
fprintf('%g\n',etime(clock,t0));

% loop through each file and create the union vector
for i=1:length(efnames)
  t_total = clock;
  
  % get the file parts
  [fpath,fname,fext] = fileparts(efnames{i});
  
  fprintf('%s:\n',fname);
  
  % set the extension and filename
  fext = sprintf('.%03i',chan);
  filename = fullfile(outdir,[fname '_union' fext]);
  
  % see if exists
  if ~exist([filename '.gz'],'file') & ~exist(filename,'file')
    % touch the file
    fid = fopen(filename,'wbl');
    fclose(fid);
    
    % get energy vector for each list
    fprintf('Calc. Power...');
    t0 = clock;
    B=single(multienergyvec(eeg{i},freqs,samplerate,width));
    fprintf('%g\n',etime(clock,t0));
    
    % calc the mean fit
    Pm = mean(double(B),2);
    
    % get the fit
    fprintf('Calc. fit...');
    t0 = clock;
    % IMPORTANT: chi_squarefit assumes that frequencies are
    % logarithmically spaced!
    all = chi_squarefit_old(freqs,Pm)';
    fprintf('%g\n',etime(clock,t0));
    
    % set the threshold
    thresh = all(:,951);
    
    % loop through each frequency and save the unionvector
    fprintf('Calc. union...');
    t0 = clock;
    unionvector = single(zeros(size(freqs,2),size(B,2)));
    
    for f = 1:size(freqs,2)
      % get the pepisode
      unionvector(f,:) = single(episodeid(B(f,:),thresh(f),3*samplerate/freqs(f),shoulder,0));
    end
    
    fprintf('%g\n',etime(clock,t0));
    
    % save the unionvector
    fid = fopen(filename,'wbl');
    fwrite(fid,size(unionvector),'int');
    fwrite(fid,unionvector,'char');
    fclose(fid);
    
    % compress it
    eval(['!gzip -f ' fullfile(outdir,[fname '_union' fext])]);
    
  end

  % give status
  fprintf('file:%d\ttime:%g\n',i,etime(clock,t_total));
end




