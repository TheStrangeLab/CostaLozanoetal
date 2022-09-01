function uvec=getp_ms(subjpath,leadno,events,durationMS,bgThreshold,datafileroot,offsetMS)
% GETU_MS - Return the Pepisode union data from a file, where times
% are in milliseconds.
%
% FUNCTION uvec=getu_ms(fpath,leadno,events,durationMS,bgThreshold,datafileroot,offsetMS)
%       subjpath = '~/eeg/free/' - the path to the directory in which
%       subjects reside
%       leadno= 34 - the electrode # to analize
%       events- events structure to analize
%       durationMS = 2000 - signal time length (from the start of
%       each event)
%       bgThreshold = 95 (optional) - threshold cutoff for the background
%       spectrum. this threshold is part of the union file name
%       (which resides in the pepisode directory) if it is other
%       than 95.
%       datafileroot = 'pepisode' (optional)- the name of the directory in
%       which the Pepisode data is located
%       offsetMS = 0 (optional)- offset at which to start the analysis
% This function returns a vector with a value of Pepisode for all
% frequencies that are in eeganalparams and for all events passed
% in. It assumes that the subject's data are in the directory with
% its subject code (e.g. BR032).

persistent lastfile

if(nargin<7) 
  offsetMS=0; % default: no offset
  if (nargin<6)
    datafileroot = 'pepisode';
    if (nargin<5)
      bgThreshold = 95;
    end
  end
end; 


% get the freqs for the size
freqs=eeganalparams('freqs');
samplerate = eegparams('samplerate');
duration = round((durationMS)*samplerate/1000);
offset = round((offsetMS)*samplerate/1000);

% allocate vector space to store all the union vectors
uvec = zeros(length(events),length(freqs));

for e = 1:length(events)

  % get the file parts
  [fpath,fname,fext] = fileparts(events(e).eegfile);

  % set the extension
  fext = sprintf('.%03i',leadno);

  % set the current file
  if (bgThreshold ==95)
  	currentfile = fullfile(subjpath,events(e).subject,datafileroot,[fname '_union' fext]);
  else
        currentfile = fullfile(subjpath,events(e).subject,datafileroot,sprintf('%s_union%i%s',fname,bgThreshold,fext));
  end

  curDir = pwd;
  cd(subjpath);

  cd(events(e).subject);

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
    here = pwd;
    begInd = findstr(lastfile,datafileroot)+length(datafileroot)+1;
    dashInd = findstr(lastfile,'_');
    endInd = dashInd(1)-1;
    lastSubj = lastfile(begInd:endInd);
    %ind = findstr('BR',lastfile);
    %lastSubj = lastfile(ind+2:ind+4);
    oldPath = [subjpath '/' lastSubj];
    cd(oldPath);
    % zip the file
    eval(['!gzip -f ' lastfile]);
    cd(here);
  end

  % open the union file
  fid = fopen(currentfile,'rb');
  if (fid <0)
     fprintf(1,'error in getp_ms --cannot find union file ');
	keyboard
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
  cd(curDir);

  % reshape to the correct size
  u = reshape(u,nrows,duration);
  % add the vector to the output vector
  uvec(e,:) = sum(u')/size(u,2);

  % set the last file to the current file
  lastfile = currentfile;
end


