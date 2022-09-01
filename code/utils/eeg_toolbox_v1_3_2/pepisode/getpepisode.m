function unionvec=getpepisode(chan,events,durationMS,offsetMS,bufferMS, varargin)
% GETPEPISODE - Return the Pepisode union data from a file. This function gives the complete
% vector of zeros and ones, you still have to do the averaging to
% come a a value for pepisode (e.g. mean(unionvec)=0.21)
%
% FUNCTION unionvec=getpepisode(chan,events,durationMS,offsetMS,
% bufferMS,varargin)
%
% INPUT ARGs:
%       leadno = 3 - the electrode #
%       events- events structure to analize
%       durationMS=2000 - signal time length in milliseconds
%       offsetMS =0 - offset at which to start in milliseconds
%       bufferMS = 0 - dummy variable to make correspondence with getphasepow
%
% OPTIONAL PARAMS:
%       'bgThreshold' = 95 - (optional) background threshold cutoff
%       for pepisode (will be used in the pepisode directory name if it
%       is other than 95)
%       'outfiledir' = 'mt' (optional), 'pepisode' is default,
%       the directory within a subject directory in which the
%       pepisode data will be written
%       'dorecompress' - when this option is included the union
%       file will be recompressed after use; the default is to
%       leave the other file open and store its filename in
%       lastfile
%       'filepath' = '~/eeg/free/' - position in which the
%       subject directory with the pepisode directory (for the
%       desired subject(s) is located
%
% OUTPUT ARGs:
%       unionvec - (Events,Freqs,Time)
%
% This function returns a vector with a value of Pepisode for all
% frequencies that are in eeganalparams. It assumes that your
% subject data is residing in the directory with the subject name
% (and this directory itself is located in the subjpath directory)

persistent lastfile

% set the defaults
bgThreshold = 95;
outfiledir = 'pepisode_standard';
dorecompress = 0;
filepath = '';

%process varargin
a = 1;
while length(varargin) > 0 & a<=length(varargin)
  switch lower(varargin{a})
   case 'bgThreshold'
    a = a+1;
    bgThreshold = varargin{a};
    a = a+1; %incremeneting a second time because the first
             %varargin is the name of this variable
             %('bgThreshold'), the second time is the argument
             %itself (e.g. 90)
   case 'outfiledir'
    a = a+1;
    outfiledir = varargin{a};
    a = a+1;
   case 'dorecompress' %don't use this option unless you have a
                       %small number of events! 
    a = a+1;
    dorecompress = 1;
   case 'filepath'
    a =a+1;
    filepath = varargin{a};
    a =a+1;
   otherwise
    error(['Error processing vararg: ' num2str(a)]);
    end
end

if bgThreshold ~= 95
  outfiledir = [outfiledir num2str(bgThreshold)];
end

% get the freqs for the size
freqs = eeganalparams('freqs');
samplerate = eegparams('samplerate');
duration = round((durationMS)*samplerate/1000);
offset = round((offsetMS)*samplerate/1000);

% allocate vector space to store all the union vectors
unionvec = zeros(length(events),length(freqs),duration);

for e = 1:length(events)
  if ~isempty(events(e).eegfile)
    % get the file parts of the eegfile, which will be used in the
    % construction of the place of the pepisode union file
    [fpath,fname,fext] = fileparts(events(e).eegfile);
    if ~isempty(filepath)
      pepisodedir = fullfile(filepath,events(e).subject,'pepisode',outfiledir); %we
      %add 'eeg.reref' to the end because the next line will remove it
    else
      pepisodedir = fullfile(fpath,'../pepisode',outfiledir);
    end
      
    % set the extension
    fext = sprintf('.%03i',chan);

    % set the current file
    currentfile = fullfile(pepisodedir,[fname '_union' fext]);
    
    % unzip the file
    if exist([currentfile '.gz'],'file')
      eval(['!gunzip ' currentfile '.gz']);
    end

    % see if is new file
    if (~isempty(lastfile) & ~strcmp(currentfile,lastfile))
      % currentfile is a new file, so compress the last file (leaving
      % files uncompressed is useful in the case a whole set of
      % events is acquired one by one for a particular subject)
      eval(['!gzip -f ' lastfile]);
      % set the lastfile variable to empty to indicate compression
      % was done
      lastfile = '';
    end

    % open the union file
    fid = fopen(currentfile,'rb','l');
    if (fid <0)
        fprintf('there is no union file\n');
	return
    end
    % get the rows and cols
    nrows = fread(fid,1,'int');
    ncols = fread(fid,1,'int');

    % read to ignore the necessary count
    thetime = offset + events(e).eegoffset;
    if thetime > 0
      fseek(fid,nrows*thetime,'cof');
    end
    % read in the amount of the data
    u = zeros(nrows,duration);
    u = fread(fid,nrows*duration,'char');

    % close the file
    fclose(fid);

    % reshape to the correct size
    if (length(u) ~= nrows*duration) %there is not enough eeg signal
      warning(sprintf('%s only %i of %i samples read for event %i, appending zeros ',events(e).eegfile,length(u),nrows*duration,e));
      u = [u; zeros(nrows*duration - length(u),1)];
    end    
    
    % reshape to the correct size
    u = reshape(u,nrows,duration);
    % add the vector to the output vector
    unionvec(e,:,:) = u;

    if dorecompress
      % compress the union file we just used
      eval(['!gzip -f ' currentfile]);
    else
      % set the last file to the current file
      lastfile = currentfile;
    end
  end  
end


