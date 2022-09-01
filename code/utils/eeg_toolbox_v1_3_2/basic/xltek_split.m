function xltek_split(eegfile,subject,outputdir,samplerate)
%XLTEK_SPLIT - Splits an ASCII XLTek file into separate channels
%
%
%
% FUNCTION:
%   xltek_split(eegfile,subject,outputdir,samplerate)
%
% INPUT ARGS:
%   eegfile = 'XLTek_session.txt';
%   subject = 'BW024';
%   outputdir = '/data/eeg/BW024/eeg.noreref';
%   samplerate = 250;
%
%




% set the defaults
if ~exist('samplerate','var')
  samplerate = 250;
end

% fix file if necessary
if ~exist([eegfile '.fixed'],'file')
  % replace AMPSAT with NaN
  fprintf('Fixing the amp. saturation...(VERY SLOW)...');
  system(['sed -e ''s/AMPSAT/NaN/g'' -e ''s/SHORT/0.0/g'' ' eegfile ' > ' eegfile '.fixed']);
  fprintf('Done\n');
end
eegfile = [eegfile '.fixed'];

% set some parameters
% ampinfo
ampinfo = [-32767 32767; -3276.7 3276.7];
ampfact = 10^6;
ampgain = calcGain(ampinfo,ampfact);
outputformat='int16';

% skip header and read in first line to get format and data/time start
fid = fopen(eegfile,'rt');
done = 0;
headerlines = 0;
while ~done
  % get a line
  s = fgetl(fid);
  
  % see if it's the EOF
  if isempty(s) | s == -1
    fprintf('Error reading file.\n')
    return
  end
  
  % see if comment
  if ~strcmp(s(1),'%')
    % is real line, so quite
    break
  end
  
  % increment header lines
  headerlines = headerlines + 1;
end
fclose(fid);

% read into cell array
txt = strread(s,'%s');

% create datevec
sdate = sscanf(txt{1},'%d/%d/%d')';
stime = sscanf(txt{2},'%d:%d:%d')';
datevec = [sdate(3) sdate(1) sdate(2) stime];

% create the output stream
filestem=fullfile(outputdir,[subject '_' datestr(datevec,'ddmmmyy_HHMM')]);

% get elec counts
numChans = length(txt)-4;  % ignore date,time,index,and word 'OFF'

% get num lines in file
[stemp,numLines] = system(['wc -l ' eegfile]);
numLines = sscanf(numLines,'%d') - headerlines;

% allocate space
eeg = zeros(numChans,numLines,'int16');
lasteeg = zeros(numChans,1);

% load in the data
fprintf('Reading in samples (%d):\n',numLines);
fid = fopen(eegfile,'r');

% read in the header lines
fprintf('Skipping %d header lines...\n',headerlines);
for i = 1:headerlines
  s = fgetl(fid);
end

for i = 1:numLines
  if i == 1 | mod(i,1000)==0
    fprintf('%d ',i);
  end
  
  ign = fscanf(fid,'%s',3);
  tempeeg = fscanf(fid,'%f',numChans);
  
  % handle saturations as just putting in the previous val
  tempeeg(isnan(tempeeg)) = lasteeg(isnan(tempeeg));
  lasteeg = tempeeg;
  
  % add it to the data
  eeg(:,i) = int16(tempeeg./ampgain);
  ign = fscanf(fid,'%s',1);
end
fclose(fid);
fprintf('%d\n',numLines);

% make sure output dir exists
if ~exist(outputdir,'dir')
  [basepath,newdir,fstem] = fileparts(outputdir);
  newdir = [newdir fstem];
 
  mkdir(basepath,newdir);
end


% loop over channels and write to files
fprintf('Writing out channel files (%d):\n',numChans);
for c = 1:numChans
  fprintf('%d ',c);
  fid = fopen([filestem '.' num2str(c,'%03d')],'wb','l');
  fwrite(fid,eeg(c,:),outputformat);
  fclose(fid);
end
fprintf('\n');


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
   
% remove the fixed file







% eeg = cell(numChans,1);
% for c = 1:numChans
%   eeg{c} = zeros(numLines,1,'single');
% end

% % make an input format for the file
% formatstr = '%*s%*s%*s';
% for i = 1:numChans
%   formatstr = [formatstr '%f'];
% end
% formatstr = [formatstr '%*s'];

% % read in the eegdata
% fprintf('Reading in samples (This could be take awhile)...\n',numLines);
% [eeg{1:numChans}] = textread(eegfile,formatstr);







% % set starting point
% fprintf('Reading in samples (%d):\n',numLines);
% ind = 1;
% while ischar(s)
%   % give status
%   if ind == 1 | mod(ind,100)==0
%     fprintf('%d ',ind);
%   end
  
%   % replace AMPSAT with NaN
%   s = strrep(s,'AMPSAT','NaN');
  
%   % split into cell array
%   txt = strread(s,'%s');
  
%   % get the eeg data out of that
%   eeg(:,ind) = str2double(txt(4:end-1));
  
%   % get the next line
%   s = fgetl(fid);
%   ind = ind + 1;
% end
% ind = ind-1;
% fprintf('%d...Done',ind);

