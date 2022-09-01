function spliteeg(eegfile, numchan, filestem, refchans, samplerate, ...
		  bandstop, ampinfo, ampfact, indataformat, outdataformat)
%SPLITEEG - Splits a scalp EEG file into separate channels
%
% Splits an EEG file into separate channels in the specified
% directory and provides the gain to convert from Raw AD values 
% to Volts (uV).
%
% At the same time you can perform a bandstop filter and
% rereference the data to an average of the specified electrodes. 
%
% When rereferencing, do not include channels recorded in bipolar
% mode.  The spliteeg function will not apply the reference to
% electrodes not in the refchans array.
%
% FUNCTION:
%   spliteeg(eegfile, numchan, filestem, refchans, samplerate, ...
%            bandstop, ampinfo, ampfact, indataformat, outdataformat)
%
% INPUT ARGUMENTS:
%   eegfile = 'dat/p300_014.dat';        % Combined eeg file
%   numchan = 64;                        % Total channels in file
%   filestem = 'dat/BR014_02Sep03_1408'; % Root name of channel
%                                        %   files
%   refchans = [1:64];                   % Channels to average for
%                                        %   reference ([] for
%                                        %   none)
%   samplerate = 256;                    % Data sample rate (default)
%   bandstop = [58 62];                  % Band stop filter range, 
%                                            use [] for no filter.
%   ampinfo = [0 4096; -5 5];            % Conversion info from raw
%                                        %   to voltage (default)
%   ampfact = 10000;                     % Amplification factor to
%                                        %   correct for (default)
%   indataformat = 'uint16';             % Format read from file
%   outdataformat = 'int16';             % Format to save to file
%

% check input args
if nargin < 10
  outdataformat = 'int16';
  
  if nargin < 9
    indataformat = 'uint16';
    
    if nargin < 8
      ampfact = 10000;
      
      if nargin < 7
        ampinfo = [0 4096; -5 5];
        
        if nargin < 6
          bandstop = [59.5 60.5];
          
          if nargin < 5
            samplerate = 256;
            
            if nargin < 4
              refchans = [];
            end
          end
        end
      end
    end
  end
end

% see if must perform amp calcs
if isempty(ampinfo)
  gain = 1;
  centerShift = 0;
else
  % perform amp conversions
  arange = abs(diff(ampinfo(1,:)));
  drange = abs(diff(ampinfo(2,:)));
  
  centerShift = arange/2 + min(ampinfo(1,:));
  ampinfo(1,:) = ampinfo(1,:) - centerShift;
  gain = calcGain(ampinfo,ampfact);
end

% determine the input databytes
if strfound(indataformat,'short') | strfound(indataformat,'16')
  databytes = 2;
elseif strcmp(indataformat,'int') | strcmp(indataformat,'uint') | strfound(indataformat,'single') | strfound(indataformat,'32')
  databytes = 4;
else
  databytes = 8;
end

% load the file 
%fprintf(['Loading ' eegfile '...']);
eeg_fid = fopen(eegfile,'r','l');

% calc the number of samples
fseek(eeg_fid,0,'eof');
totalsamps = ftell(eeg_fid)/databytes;
fseek(eeg_fid,0,'bof');
numsamps = totalsamps/(numchan);

% loop to load from file so it takes less memory
stepsize = 1000000;
fprintf('Loading %g million samples(%d)...\n',totalsamps/stepsize,totalsamps);

totalread = 0;
eeg = int16(zeros(totalsamps,1));
while totalread < totalsamps
  sampsleft = totalsamps - totalread;
  if sampsleft < stepsize
    toread = sampsleft;
  else
    toread = stepsize;
  end
  fprintf('%g ',totalread/stepsize);
  eeg(totalread+1:totalread+toread) = int16(fread(eeg_fid,toread,indataformat)./gain);
  totalread = totalread + toread;
end
fprintf('%d...Done\n',totalread);

% close the main file
fclose(eeg_fid);

%eeg = int16(fread(fid,indataformat));
%fclose(fid);

% resize the data
eeg = reshape(eeg,numchan,length(eeg)/numchan);

% see if rerference
if length(refchans) > 0
  fprintf('\nRereferencing...');
  
  eeg_avg = round(mean(eeg(refchans,:),1));
  
  % save the avg reference
  % make the filename
  chanfile = sprintf('%s.avg', filestem);
  
  % open and write the file
  fid = fopen(chanfile,'w','l');
  fwrite(fid,eeg_avg,outdataformat);
  fclose(fid);
end

fprintf('\nProcessing %d channels:\n', numchan);
for i = 1:numchan
  fprintf('%i ', i);
  
  % make the filename
  chanfile = sprintf('%s.%03i', filestem, i);
  
  % set the output
  eout = double(eeg(i,:));

  % see if rereference
  if ismember(i,refchans)
    % subtract the reference, which centers at zero
    eout = eout - eeg_avg;
  else
    % center around zero
    eout = eout - centerShift;
  end
  
  % see if filter
  if length(bandstop) > 0
    % filter it
    eout = buttfilt(eout,bandstop,samplerate,'stop',1);
  end
  
  % open and write the file
  fid = fopen(chanfile,'w','l');
  fwrite(fid,eout,outdataformat);
  fclose(fid);
end

% write out params.txt file
[pathstr,name,ext] = fileparts(filestem);
paramfile = fullfile(pathstr,'params.txt');
fid = fopen(paramfile,'w','l');
fprintf(fid,'samplerate %d\ndataformat ''%s''\ngain %g\n',samplerate,outdataformat,gain);
fclose(fid);

fprintf('\n\n');

