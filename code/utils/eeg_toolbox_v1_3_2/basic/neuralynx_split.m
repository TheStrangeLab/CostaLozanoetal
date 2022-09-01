function neuralynx_split(sessID,subject,outDir)
%Neuralynx_split - Splits a Neuralynx .eeg file into separate channels
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
%   neuralynx_spliteeg(ncsdir,subject,outputdir)
%
% INPUT ARGUMENTS:
%    sessID = 'UP011_18Jul07_1149neuralynx'
%    subject = 'UP011'
%    outDir = 'eeg.noreref.neuralynx'
%

% Define directories
basedir=fullfile('/data/eeg',subject);
ncsdir=fullfile(basedir,'raw',sessID);
outputdir=fullfile(basedir,outDir);

% Define outputfilename
if strcmp(sessID(end-8:end),'neuralynx')
   outputfilename=sessID(1:end-9);
else
   outputfilename=sessID;
end


outdataformat='int16';

eegSampleRate=2000;
resampledDT=1e6/eegSampleRate;

if ~exist(outputdir,'dir')
  mkdir(outputdir)
end

paramsfile=fullfile(outputdir,'params.txt');
if ~exist(paramsfile,'file')
  fid=fopen(paramsfile,'w');
  fprintf(fid,'samplerate %d\n',eegSampleRate);
  fprintf(fid,'dataformat ''single''\n');
  fclose(fid);
end


for n=1:96
  possibleFilenames{n}=sprintf('CSC%d',n);
end
possibleFilenames={possibleFilenames{:},'RefA1','RefA2','RefA3','RefB1', ...
                   'RefB2','RefB3','RefC1','RefC2','RefC3','RefD1','RefD2', ...
                   'RefD3'};

curChan=1;
jackFlag=0;

for filenum=1:length(possibleFilenames)

  fName=fullfile(ncsdir,[possibleFilenames{filenum} '.Ncs']);
  if ~exist(fName),  fName=fullfile(ncsdir,[possibleFilenames{filenum} '.ncs']);end
  if ~exist(fName), continue;
  end

  fprintf('Loading %s\n',fName);
  [data,info,allTS]=load_ncs(fName);
  
   % outName=fullfile(outputdir,sprintf('%s_%s.%.03i',subject, ...
   %                                  datestr(datenum(info.fileStartTime), ...
   %                                          'ddmmmyy_HHMM'),curChan));
   % Change this to call output files same names as inputfiles
   outName=fullfile(outputdir,[outputfilename sprintf('.%.03i',curChan)]);

  if jackFlag==0
    % jackFile=fullfile(outputdir,sprintf('%s_%s.jacksheet.txt',subject, ...
    %                                  datestr(datenum(info.fileStartTime), ...
    %                                          'ddmmmyy_HHMM')));
    jackFile=fullfile(outputdir,[outputfilename '.jacksheet.txt']);

    jackFid=fopen(jackFile,'w');
    jackFlag=1;
  end

    dataTimes=linspace(info.firstSampleTime,info.lastSampleTime,length(data));  


  resampledTimes=info.firstSampleTime:resampledDT:info.lastSampleTime;
  dataVolts=-double(data).*info.bitVolts;
  
  clear data;pack


  numChunks=10;
chunkBoundaries=round(linspace(1,length(resampledTimes),numChunks+1));
  %here we are splitting into chunks so that we don't run out of memory in interp1 on large datasets. i know this is so fucking stupid but kareem and i couldn't think of anything better
lfp=zeros(size(resampledTimes));
  for i=1:numChunks

    resampIdx=chunkBoundaries(i):chunkBoundaries(i+1);
    dataStartIdx=find(dataTimes<=resampledTimes(chunkBoundaries(i)),1,'last');
    dataEndIdx=find(dataTimes>=resampledTimes(chunkBoundaries(i+1)),1,'first');

    dataTimesTmp=dataTimes([dataStartIdx:dataEndIdx]);
    dataVoltsTmp=dataVolts([dataStartIdx:dataEndIdx]);
    resampTimesTmp=resampledTimes(resampIdx);
    lfp(resampIdx)=interp1(dataTimesTmp,dataVoltsTmp,resampTimesTmp,'linear');

  end


  pack

  fprintf('\t writing to %s\n',outName);
  outfid=fopen(outName,'w','l');
  c=fwrite(outfid,lfp,'single');
  fclose(outfid);

  fprintf(jackFid,'%d\t%s\n',curChan,possibleFilenames{filenum});
  curChan=curChan+1;


end


fclose(jackFid);

% Now create pulses.sync.txt file that records samples of all pulses

eventsFile=fullfile(ncsdir,'Events.nev');
eventsTS=load_nev(eventsFile);

% Normalize events by initial time of recording
eventsTS=eventsTS-info.firstSampleTime;
eventsSamples=round(eventsTS/resampledDT)+1;

% Write the pulse samples
% pulsefile=fullfile(outputdir,sprintf('%s_%s.sync.txt',subject, ...
%                                      datestr(datenum(info.fileStartTime), ...
%                                              'ddmmmyy_HHMM')));
pulsefile=fullfile(outputdir,[outputfilename '.sync.txt']);


pulseFid=fopen(pulsefile,'w','l');
fprintf(pulseFid,'%d\n',eventsSamples);
fclose(pulseFid);




function [upstrokes,downstrokes]=load_nev(filename)

% [uptimes, downtimes]=load_nev(filename)
% This function reads in the specified Neuralynx .Nev file, and outputs the
% upstrokes and downstrokes from it.
%
% From our setup at UCLA, we found that the TTL value of -2 corresponds to an
% upstroke, and -4 corresponds to a downstroke.  This function makes that assumption!

%revised 8/12: sometimes upstrokes seem to be -1, and downs seem to be
%-3. changed code accrdingly

% revised 7/23 for Penn neuralynx system - no upstrokes seem to be 1 and
% downstrokes seem to be 0.  there are occasional ttlvalues of 3, but we'll
% ignore these as they don't seem to be related to the actual sync pulses

fid = fopen(filename,'r','l');

fseek(fid,16*1024,'bof');

upstrokes=[];
downstrokes=[];
ttlValues=[];

while ~feof(fid)
  pktStart=fread(fid,1,'int16');
  
  pktId=fread(fid,1,'int16');
  pktDataSize=fread(fid,1,'int16');
  timeStamp=fread(fid,1,'int64');
  eventId=fread(fid,1,'int16');
  ttlValue=fread(fid,1,'int16');
  ttlValues=[ttlValues ttlValue];
  crc=fread(fid,1,'int16');
  dummy=fread(fid,1,'int32');
  extra=fread(fid,8,'int32');
  eventString=char(fread(fid,128,'char')');

  if isempty(ttlValue)
    disp(sprintf('%s ended prematurely!',filename));
    break
  end  
  
%   if ttlValue== -2 || ttlValue==-1
%     upstrokes=[upstrokes timeStamp];
%   elseif ttlValue==-4 || ttlValue==-3
%     downstrokes=[downstrokes timeStamp];
%   end
  if ttlValue== -2 || ttlValue==-1 || ttlValue==1
    upstrokes=[upstrokes timeStamp];
  elseif ttlValue==-4 || ttlValue==-3 || ttlValue==0
    downstrokes=[downstrokes timeStamp];
  end
  
end

fclose(fid);







