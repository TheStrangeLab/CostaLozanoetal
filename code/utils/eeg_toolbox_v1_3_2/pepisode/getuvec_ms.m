function unionvec=getuvec_ms(subjpath,leadno,events,durationMS,bgThreshold,outfiledir,offsetMS,freqs)
% GETUVEC_MS - Return the Pepisode union data from a file. This function gives the complete
% vector of zeros and ones, you still have to do the averaging
%
% FUNCTION unionvec=getuvec_ms(subjpath,leadno,events,durationMS,bgThreshold,outfiledir,offsetMS)
%       subjpath = '~/eeg/scalp_free/' - home directory of the experiment in which the subject directories reside
%       leadno = 3 - the electrode #
%       events- events structure to analize
%       durationMS=2000 - signal time length in milliseconds
%       bgThreshold = 95 - (optional) background threshold cutoff for pepisode
%       outfiledir = 'pepisode_standard' (optional), 'pepisode' is default, the directory within a subject directory in which the pepisode data will be written
%       offsetMS =0 (optional) - offset at which to start in
%       milliseconds
%       freqs =(2^(1/8)).^(8:56) (optional) - set of freqs for
%       which Pepisode should be computed 
% This function returns a vector with a value of Pepisode for all
% frequencies that are in eeganalparams. It assumes that your
% subject data is residing in the directory with the subject name

persistent lastfile

if (nargin<8)
  freqs = eeganalparams('freqs');
  if(nargin<7) 
    offsetMS=0;
    if (nargin < 6)
     outfiledir = 'pepisode_standard';
     if (nargin < 5)
        bgThreshold = 95;
     end
    end 
  end; % default: no offset
end

% get the freqs for the size
samplerate = eegparams('samplerate');
duration = round((durationMS)*samplerate/1000);
offset = round((offsetMS)*samplerate/1000);

% allocate vector space to store all the union vectors
unionvec = zeros(length(events),length(freqs),duration);

for e = 1:length(events)
  if ~isempty(events(e).eegfile)
    % get the file parts
    [fpath,fname,fext] = fileparts(events(e).eegfile);

    % set the extension
    fext = sprintf('.%03i',leadno);

    % set the current file
    if (bgThreshold ==95)
  	currentfile = fullfile(subjpath,events(e).subject,'pepisode',outfiledir,[fname '_union' fext]);
    else
        currentfile = fullfile(subjpath,events(e).subject,'pepisode',sprintf('%s%i',outfiledir,bgThreshold),sprintf('%s_union%s',fname,fext));
    end


    % unzip the file
    if exist([currentfile '.gz'],'file')
      eval(['!gunzip ' currentfile '.gz']);
    end

    % see if is new file
    if (~isempty(lastfile) & ~strcmp(currentfile,lastfile))
      % is new file, so compress the last file; we have to deal with
      % the possibility that the last file was from another subject;
      % this assumes that every subject EEG file starts with the
      % subject's code
      %here = pwd;
      %begInd = findstr(lastfile,outfiledir)+length(outfiledir)+1;
      %dashInd = findstr(lastfile,'_');
      %endInd = dashInd(1)-1;
      %lastSubj = lastfile(begInd:endInd);
      %ind = findstr('BR',lastfile);
      %lastSubj = lastfile(ind+2:ind+4);
      %oldPath = [subjpath '/' lastSubj];
      %cd(oldPath);
      % zip the file
      eval(['!gzip -f ' lastfile]);
    end

    % open the union file
    fid = fopen(currentfile,'rb');
    if (fid <0)
        warning('there is no union file\n')
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
    %cd(curDir);

    % reshape to the correct size
    if (length(u) ~= nrows*duration) %there is not enough eeg signal
      warning(sprintf('%s only %i of %i samples read for event %i, appending zeros ',events(e).eegfile,length(u),nrows*duration,e));
      u = [u; zeros(nrows*duration - length(u),1)];
    end
    u = reshape(u,nrows,duration);
    % add the vector to the output vector
    unionvec(e,:,:) = u;

    % set the last file to the current file
    lastfile = currentfile;
  end %if there is an eegfile
end


