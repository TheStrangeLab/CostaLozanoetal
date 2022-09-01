function nicolet_split(eegfile,subject,outputdir)
%nicolet_split - Splits a Nicolet .eeg file into separate channels
%
% Splits an EEG file into separate channels in the specified
% directory.  Each EEG file must have an accompanying .bni file
% containing header information Writes a params.txt and jacksheet.txt
% file based on the parameters from the .bni file.
%
% Output file base is created as a function of the EEG recording start
% time.
%
%
% FUNCTION:
%   nicolet_split(eegfile,subject,output_dir)
%
% INPUT ARGUMENTS:
%   eegfile = 'data.001';        % Combined eeg file
%   subject = 'UP001'
%   output_dir='/data/eeg/UP001/eeg.noreref'
%

outdataformat='int16';

% set the eeg and bni files
bnifile=[eegfile '.bni'];
if ~exist(eegfile,'file')
  eegfile = [eegfile '.eeg'];
  if ~exist(eegfile,'file')
    fprintf('EEGfile not found (with or without .eeg).\n')
    return
  end
end

% check whether the output directory exists, and if not, create it
if ~exist(outputdir,'dir')
  mkdir(outputdir);
end

disp(['Loading '  bnifile '...']);
hdr=EEGbni(bnifile,eegfile);

disp(['Loading '  eegfile '...']);
eeg=nicoletReadEEG(eegfile,hdr);

%filestem=fullfile(outputdir,[subject '_' datestr(hdr.datevec,'ddmmyy_HHMMSS')]);
filestem=fullfile(outputdir,[subject '_' datestr(hdr.datevec,'ddmmmyy_HHMM')]);

fprintf('Processing %d channels:\n', hdr.nchan);
curChan = 0;

for i = 1:hdr.maxchan
  if isempty(hdr.labels{i})
    continue
  end

  curChan = curChan + 1;

  fprintf('%i ', i);
  
  % make the filename
  chanfile = sprintf('%s.%03i', filestem, i);
  
  % open and write the file
  fid = fopen(chanfile,'w','l');
  fwrite(fid,eeg(curChan,:),outdataformat);
  fclose(fid);
end
fprintf('\n');

% write out params.txt file
[pathstr,name,ext] = fileparts(filestem);
paramfile = fullfile(pathstr,'params.txt');
fid = fopen(paramfile,'w','l');
fprintf(fid,'samplerate %.2f\ndataformat ''%s''\ngain %g\n',hdr.rate,outdataformat,hdr.sens);
fclose(fid);

% new params matching base name
paramfile = [filestem '.params.txt'];
fid = fopen(paramfile,'w');
fprintf(fid,'samplerate %d\ndataformat ''%s''\ngain %g\n',hdr.rate,outdataformat,hdr.sens);
fclose(fid);


%print the labels for each channel
fopen(fullfile(outputdir,'jacksheet.txt'),'w','l');
for i=1:size(hdr.labels,1)
  fprintf(fid,'%d %s\n',i,hdr.labels{i});
end
fclose(fid);


function raw = nicoletReadEEG(inputfile,hdr)
%reads all raw data from nicolet .eeg file
%this code is copied from the Litt Lab's toolbox
%JJ modified it for use with our toolbox.

%check to see if this the old Nicolet compressed file format
if strcmp(hdr.hdrtype,'TAG')
  error('JJ: no support for the old Nicolet compressed format for now')
end

totalsamps = hdr.npts * hdr.nchan;

fprintf('Loading %d samples...\n',totalsamps);

fid = fopen(inputfile,'r','l');

% loop to load from file so it takes less memory
stepsize = 1000000;
totalread = 0;
raw = int16(zeros(totalsamps,1));
while totalread < totalsamps
  sampsleft = totalsamps - totalread;
  if sampsleft < stepsize
    toread = sampsleft;
  else
    toread = stepsize;
  end
  fprintf('%d ',totalread);
  raw(totalread+1:totalread+toread) = int16(fread(fid,toread,'int16'));
  totalread = totalread + toread;
end
fprintf('%d...Done\n',totalread);

% reshape for processing
raw = reshape(raw,hdr.nchan, hdr.npts);

% close the main file
fclose(fid);

%kx = fread(fid,[hdr.nchan hdr.npts],'int16');



function hdr = EEGbni(bnifile,datafile)
%EEGbni - gets parameters from nicolet .bni header file
%   HdrStruct = EEGbni(filename)
%       HdrStruct - and eeg header structure (see EEGHdrStruct)

global DEBUGdisp

hdr = struct('npts',[],'rate',[],'nchan',[],'sens',[],'datevec',[],...
    'labels',[],'hdrtype',[],'nextfile',[],'events',[],'mpginfo',[],'hdrbytes',[],'etime',[],'montage',[]);
sdate = [];
stime = [];
fid = fopen(bnifile, 'r','l');
if fid == -1, error('File not found'); end
while 1,
   txt = fgetl(fid);%disp(txt);
   if ~isstr(txt), break; end
   if ~isempty(findstr('Filename',txt)), fname = sscanfname(txt,'Filename = %s'); end
   if ~isempty(findstr('Date = ',txt)), sdate = sscanf(txt,'Date = %d/%d/%d')'; end
   if ~isempty(findstr('Time =',txt)), stime = sscanf(txt,'Time = %d:%d:%d')'; end
   if ~isempty(findstr('Rate',txt)), hdr.rate = sscanf(txt, 'Rate = %f Hz'); end
   if ~isempty(findstr('NchanFile',txt)), hdr.nchan = sscanf(txt, 'NchanFile = %d'); end
   if (~isempty(findstr('UvPerBit',txt))) & (isempty(hdr.sens)), 
      hdr.sens = sscanf(txt, 'UvPerBit = %f'); 
   end
   %josh changed this next line from:
   %   if ~isempty(findstr('MontageRaw',txt)) & (~exist('labels')),
   % to this:
   % pbs changed MontageRaw to MontageGaped
   %
      if ~isempty(findstr('MontageGaped =',txt))
       if DEBUGdisp, disp('MontageGaped'); end
       try,
         I = findstr(txt,',');
         hdr.maxchan = max(find(diff(I)>1))+1;
         hdr.labels = cell(hdr.maxchan,1);
         %hdr.labels = strvcat(hdr.labels,txt(14:I(1)-1));
         ind1 = 16;
	 ind2 = I(1)-1;
         if ind1 < ind2
           hdr.labels{1} = txt(ind1:ind2);
         end
         for chan = 1:hdr.maxchan-1,
             ind1 = I(chan)+1;
             ind2 = I(chan+1)-1;
             if ind1 ~= ind2
               % add the label
               hdr.labels{chan+1} = txt(ind1:ind2);
             end
             %hdr.labels = strvcat(hdr.labels,txt(ind1:ind2));
         end
         
      catch, disp(lasterr);
      end
   end
   if ~isempty(findstr('[Events]',txt)),
       if DEBUGdisp, disp('[Events]'); end
       nmpg = 0;
       txt = fgetl(fid);
       while isempty(findstr('NextFile',txt)) & isstr(txt),
           % find montage settings
           if ~isempty(findstr('Montage:',txt)),
               if DEBUGdisp, disp('Montage'); end
               I = findstr('Montage:',txt);
               txt = txt(I(1):end);
               I = findstr(txt,'Selected');
               if ~isempty(I), 
                   if DEBUGdisp, disp('Montage: Selected Lines:'); end
                   mtgfmt = 'Montage: Selected Lines: %d';
                   nlines = sscanf(txt,mtgfmt);
               else, 
                   if DEBUGdisp, disp('Montage: Lines:'); end
                   mtgfmt = 'Montage: %s Lines: %d';
                   A = sscanf(txt,mtgfmt);
                   nlines = A(end);
               end
               for k = 1:nlines,
                   if DEBUGdisp, disp('Lines:'); end
                   txt = fgetl(fid);
                   hdr.montage = strvcat(hdr.montage,txt); %can parse this later for clipping
               end
           else,
               if DEBUGdisp, disp('Events'); end
               hdr.events = strvcat(hdr.events,txt);
           end
           txt = fgetl(fid);
       end
       if ~isempty(findstr('NextFile',txt)), hdr.nextfile = sscanf(txt,'NextFile = %s'); end
   end
end

fclose(fid);
hdr.datevec = [sdate(3) sdate(1) sdate(2) stime];
%get number of data points
if nargin < 2, [pa,fname,ext] = fileparts(fname); datafile = [fname ext]; end

D = dir(deblank(datafile));
hdr.npts = D.bytes;
hdr.npts = hdr.npts./(2*hdr.nchan);
hdr.npts = floor(hdr.npts);
hdr.etime = hdr.npts/hdr.rate;


function filestr = sscanfname(filestr,leadstr)

%SSCANFNAME - sscanf function to deal with filenames that include spaces
%   filepath = sscanfname(filestr,leadstr)
%   filestr - string to be parsed
%   leadstr - optional string to be separated
%   example:
%   fname = sscanfname('Filename = c:\dir name\file name.ext','Filename = ')
%   Will return fpath = 'c:\dir name\file name.ext'
%   Normally sscanf cannont handle blanks in the middle of the string so
%        fpath = sscanf('Filename = c:\dir name\file name.ext','Filename = %s')
%   will return the string fpath = 'c:\dir' only
%
%   created: 2/7/03 (sdc)

if nargin < 2, leadstr = []; end

%first replace blanks in strings with dummy char
filestr = strrep(filestr,' ','$');
leadstr = strrep(leadstr,' ','$');
%now parse string with sscanf
leadstr = [leadstr '%s'];
filestr = sscanf(filestr, leadstr);
%now go back and insert blanks
filestr = strrep(filestr, '$',' ');
