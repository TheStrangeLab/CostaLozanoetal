function spikevec = get_spikevec_ms(unitno,events,durationMS,offsetMS,outfileDir)
% GET_SPIKEVEC_MS - Return the vector with zeros and ones for whether ...
%       a spike has been found at least sample in time
%
% FUNCTION spikevec =
% get_spikevec_ms(unitno,events,durationMS,offsetMS,outfileDir)
% 
% INPUT ARGs
%  unitno = 2               % number of the cell you want to look at
%  events = events(8:12)    % events struct to extract
%  durationMS = 1000        % duration in ms of the window you want
%                           % to look at
%  offsetMS = 0             % offset at which to start in milliseconds
%  outfileDir = 'spikes'    % directory in which to look for the
%                           % cell files
%
% OUTPUT Args
%  spikevec(trials,time) - matrix with 0 when there is no spike at
%  that sample, 1 when there is a spike
%
%

if ~exist('outfileDir','var')
  outfileDir = 'spikes';
end

fulloutfileDir = fullfile('/data/eeg/',events(1).subject,outfileDir);
samplerate = eegparams('samplerate',fulloutfileDir);
duration = round((durationMS)*samplerate/1000);
offset = round((durationMS)*samplerate/1000);
%allocate vector space to store the union vectors
spikevec = zeros(length(events),duration);
lastfile = '';
for e = 1:length(events)
  % open the relevant cell
  datFile = fullfile('/data/eeg/',events(e).subject,sprintf('unit_%s.%03d',events(e).unitFile,unitno));
  if ~strcmp(datFile,lastfile)
    % zip the previous file     
    eval(['!gzip -f ' lastfile]);
    % we have to open and unzip the cell
    % if necessary, unzip the file
    if exist([datFile '.gz'],'file')
      eval(['!gunzip ' datFile '.gz']);
    end
  end
  fid = fopen(datFile,'rb');
  if (fid <0)
    warning('there is no union file\n')
    return
  end
    
  % read to ignore the necessary count
  thetime = offset + events(e).unitoffset;
  if thetime > 0
    fseek(fid,thetime,'cof');
  end    
  % read in the amount of the data
  spikevec(e,:) = fread(fid,duration,'char');
  fclose(fid);
  
  
  % set the last file to the current file
  lastfile = currentfile;
    


end
