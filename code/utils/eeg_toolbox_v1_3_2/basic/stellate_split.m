function basename = stellate_split(eegfile,subject,outputdir)
% STELLATE_SPLIT - Process a continuous file exported in ASCII (text)
% format from Stellate (used at UCLA).
%
% FUNCTION:
%    stellate_split(eegfile,subject,outputdir)
%
% INPUT ARGS:
%  eegfile = 'U397_27Mar1940.txt' (required)
%  subject = 'U397' (required)
%  outputdir = '/data/eeg/U397/eeg' - directory to put split out
%  channel files (default: '.')
%
% OUTPUT ARGS:
%  basename - the basename determined from the subject and file
%

% check input vars
if ~exist('subject','var')
  error('You must supply a subject ID.');
end

if ~exist('outputdir','var')
  outputdir = '.';
end

%check whether outputdir exists, if not, create it
if ~exist(outputdir,'dir')
  mkdir(outputdir);
end

% months
month = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

outputformat = 'single';
ampgain = 1;

% open the file
fid = fopen(eegfile,'r','l');
if fid == -1
  error('eegfile %s not found',eegfile);
end

% read the header to get info for basename
while 1,
  txt = fgetl(fid);
  if ~isempty(findstr('"TIME"',txt)), break, end
  %if ~isempty(findstr('"Patient name"',txt)), hdr.patNum = sscanf(txt,'"Patient name", "%s"'); hdr.patNum = strrep(hdr.patNum,'"',''); end
  if ~isempty(findstr('"Date recorded"',txt)), hdr.date = sscanf(txt,'"Date recorded", %d-%d-%d')'; end
  if ~isempty(findstr('"Time recorded"',txt)), hdr.time = sscanf(txt,'"Time recorded", %s:%s:%s'); end
  if ~isempty(findstr('uV',txt)), [chanNumLoc,eegtype,units,rate] = strread(txt,'"#%s%s%s%s"','delimiter',','); hdr.rate = str2num(strrep(cell2mat(rate),' Hz"','')); [nchan,loc] = strread(chanNumLoc{:},'%d%s','delimiter',' '); hdr.nchan = nchan; end
end

% construct basename
day = num2str(hdr.date(3));
themonth = hdr.date(2);
year = num2str(hdr.date(1));
time = [char(hdr.time(1:2)) char(hdr.time(4:5))];
basename = [subject '_' day month{themonth} year(end-1:end) '_' time];

% Give them some file info
fprintf('EEG File Information:\n')
fprintf('---------------------\n')
fprintf('Sample Rate = %d\n', hdr.rate);
fprintf('Start of recording = %d/%d/%d: %s\n',str2num(day),themonth,str2num(year),hdr.time);
fprintf('Number of channels = %d\n', hdr.nchan);
fprintf('Basename = %s\n', basename);
fprintf('\n');

% account for the timepoint column being split into 3 columns
nchanplus = hdr.nchan + 3;

% load the rest of the file
fprintf('Processing %d channels...\n',hdr.nchan);
% read the data in steps
fprintf('Loading data in chunks: ');
stepsize = 100000; % lines
dat = {0};
datsize = -1;
while 1,
  datsize = datsize + size(dat{1},1);
  fprintf('%d ',datsize);
  dat = textscan(fid,['%d:%d:%.3f' repmat('%.6f',1,hdr.nchan)],stepsize,'Delimiter',',');
  % write/append the data to file for every channel
  for c = 4:nchanplus
    fidw = fopen(fullfile(outputdir,[basename '.' num2str(c-3,'%03d')]),'ab','l');
    fwrite(fidw,dat{c}(:),outputformat);
    fclose(fidw);
  end
  % break if at the end of the file
  if size(dat{1},1) ~= stepsize
    break
  end
end
datsize = datsize + size(dat{1},1);
fprintf('%d\n',datsize);

fclose(fid);

% write out params.txt file
paramfile = fullfile(outputdir,'params.txt');
fid = fopen(paramfile,'w');
fprintf(fid,'samplerate %d\ndataformat ''%s''\ngain %g\n',hdr.rate,outputformat,ampgain);
fclose(fid);
           
% write out new params.txt file
paramfile = fullfile(outputdir,[basename '.params.txt']);
fid = fopen(paramfile,'w');
fprintf(fid,'samplerate %d\ndataformat ''%s''\ngain %g\n',hdr.rate,outputformat,ampgain);
fclose(fid);
