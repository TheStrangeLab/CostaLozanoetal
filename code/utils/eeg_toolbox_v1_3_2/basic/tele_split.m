function tele_split(fname,filestem,samplerate)
%TELE_SPLIT - Split out Telefactor ASCII EEG data.
%
%
% FUNCTION:
%   tele_split(eegfile,filestem,samplerate)
%
% INPUT ARGS:
%   eegfile = 'CP001_01Nov05_1700.txt';   % The eeg file to read
%   filestem = '../eeg.noreref/CP001_01Nov05_1700' % Where to put
%                                 % the data and the base file stem
%   samplerate = 200;  % The sample rate of the data
%

if ~exist('samplerate','var')
  samplerate = 200;
end

outputdir = fileparts(filestem);

% set some parameters
% ampinfo
ampinfo = [-32767 32767; -3276.7 3276.7];
ampfact = 10^6;
ampgain = calcGain(ampinfo,ampfact);
outputformat='int16';

maxbuf=10000; % max. number of timepoints to read in at a time
in=fopen(fname,'r');

buf=fscanf(in,'%i',2); % read in the header
nChan=buf(2);

% open all the output files
for c=1:nChan
  outfname=sprintf('%s.%03i',filestem,c);
  out(c)=fopen(outfname,'w','l');
  if(out(c)==-1)
    fprintf(1,'\nError: Can''t open file %s\n',outfname);
    return;
  end
end % for through channels

okay=1;
fprintf('Extracting data...');
while(okay)
  fprintf('.');

  buf=fscanf(in,'%g',[nChan,maxbuf]);
 
  if(isempty(buf)) 
    okay=0;
  % write to each channel:
  else  
    % apply the gain and write it out
    buf = int16(buf./ampgain);
    for c=1:nChan
      fwrite(out(c),buf(c,:),outputformat); 
    end; 
  end;
end % looping through the file

fprintf('Done\n')

for c=1:nChan
  fclose(out(c)); 
end % close all output files

fclose(in);

% write out params.txt file
paramfile = fullfile(outputdir,'params.txt');
fid = fopen(paramfile,'w');
fprintf(fid,'samplerate %d\ndataformat ''%s''\ngain %g\n',samplerate,outputformat,ampgain);
fclose(fid);

% new params matching base name
paramfile = [filestem '.params.txt'];
fid = fopen(paramfile,'w');
fprintf(fid,'samplerate %d\ndataformat ''%s''\ngain %g\n',samplerate,outputformat,ampgain);
fclose(fid);
   


