function split_ncs(sessID,subject,secondSess)
%split_ncs - Splits a Neuralynx .eeg file that contains a large 
%            break into two separate .eeg files, one with the 
%            original name, and one with the new name, secondSess
%
%
% FUNCTION:
%   split_ncs(sessID,subject,secondSess)
%
% INPUT ARGUMENTS:
%    sessID = 'UP011_19Jul07_1537'
%    subject = 'UP011'
%    secondSess = 'UP011_19Jul07_1540'
%

% Define directories
basedir=fullfile('/data/eeg',subject);
ncsdir=fullfile(basedir,'raw',sessID);
secondncsdir=fullfile(basedir,'raw',secondSess);

% Make this second directory
if ~exist(secondncsdir,'dir')
     mkdir(secondncsdir)
end

% Define bytes and other constants
bytesPerRecord=1024+20;
headerBytes=16*1024;
samplesPerRecord=512;

sampRate=32556;



% Set threshold for break time - timestamps that are greater than this much 
% time apart indicate where we should split the file
brkThrsh=1e6;   % One second

for n=1:96
  possibleFilenames{n}=sprintf('CSC%d',n);
end
possibleFilenames={possibleFilenames{:},'RefA1','RefA2','RefA3','RefB1', ...
                   'RefB2','RefB3','RefC1','RefC2','RefC3','RefD1','RefD2', ...
                   'RefD3'};

curCh=0;

for filenum=1:length(possibleFilenames)

  fName=fullfile(ncsdir,[possibleFilenames{filenum} '.Ncs']);
  if ~exist(fName),  fName=fullfile(ncsdir,[possibleFilenames{filenum} '.ncs']);end
  if ~exist(fName), continue;
  end
  
  fNameSecond=fullfile(secondncsdir,[possibleFilenames{filenum} '.ncs']);

  % Grab stats
  info=stat_ncs(fName);

  % Grab header information for writing later and all of data
  fid=fopen(fName,'r','l');
  headerString=fread(fid,headerBytes,'*char');

  allTS=zeros(1,info.numRecords);
  data=zeros(1,info.numRecords*samplesPerRecord);

  for record=0:info.numRecords-1
    ts=fread(fid,1,'int64');

    tmp=fread(fid,3,'int32');
    curCh=tmp(1);
    sampFreq=tmp(2);
    numValidSamp=tmp(3);

    thisRecordData=fread(fid,512,'*int16');
    allTS(record+1)=ts;
    data(record*samplesPerRecord+1:(record+1)*samplesPerRecord)=thisRecordData;
  end

  fclose(fid);

  % Find point where data breaks over long period of time
  diffTS=diff(allTS);
  breakPt=find(diffTS>brkThrsh);
  dataBrk=breakPt*samplesPerRecord;

  % Now divide data and timestamps into first and second

  allTSFirst=allTS(1:breakPt);
  allTSSecond=allTS(breakPt+1:end);

  dataFirst=data(1:dataBrk);
  dataSecond=data(dataBrk+1:end);

  numRecFirst=length(allTSFirst);
  numRecSecond=length(allTSSecond);

  % Now write two halves of file into two separate files

  % First file

  fid=fopen(fName,'w','l');
  fwrite(fid,headerString,'char');

  for record=0:numRecFirst-1
    fwrite(fid,allTSFirst(record+1),'int64');

    fwrite(fid,curCh,'int32');
    fwrite(fid,sampFreq,'int32');
    fwrite(fid,numValidSamp,'int32');

    fwrite(fid,dataFirst(record*samplesPerRecord+1:(record+1)*samplesPerRecord),'int16');
  
  end

  fclose(fid);

  % Second file

  fid=fopen(fNameSecond,'w','l');
  fwrite(fid,headerString,'char');

  for record=0:numRecSecond-1
    fwrite(fid,allTSSecond(record+1),'int64');

    fwrite(fid,curCh,'int32');
    fwrite(fid,sampFreq,'int32');
    fwrite(fid,numValidSamp,'int32');

    fwrite(fid,dataSecond(record*samplesPerRecord+1:(record+1)*samplesPerRecord),'int16');
  
  end

  fclose(fid);


end


% Now read in events file and rewrite a new one
eventsFile=fullfile(ncsdir,'Events.nev');
eventsFileSecond=fullfile(secondncsdir,'Events.nev');

d=dir(eventsFile);
bytesPerEvRecord=184;

numEventsRec=(d.bytes-headerBytes)/bytesPerEvRecord;


% Read in events file
fid=fopen(eventsFile,'r','l');


eventHeader=fread(fid,headerBytes,'*char');

pktStart=zeros(1,numEventsRec);
pktId=zeros(1,numEventsRec);
pktDataSize=zeros(1,numEventsRec);
timeStamp=zeros(1,numEventsRec);
eventId=zeros(1,numEventsRec);
ttlValues=zeros(1,numEventsRec);
crc=zeros(1,numEventsRec);
dummy=zeros(1,numEventsRec);
extra=zeros(1,numEventsRec);
eventString=cell(1,numEventsRec);

for record=0:numEventsRec-1

  pktStart(record+1)=fread(fid,1,'int16');
  pktId(record+1)=fread(fid,1,'int16');
  pktDataSize(record+1)=fread(fid,1,'int16');
  timeStamp(record+1)=fread(fid,1,'int64');
  eventId(record+1)=fread(fid,1,'int16');
  ttlValues(record+1)=fread(fid,1,'int16');
  crc(record+1)=fread(fid,1,'int16');
  dummy(record+1)=fread(fid,1,'int32');
  extra(record*8+1:(record+1)*8)=fread(fid,8,'int32');
  eventString{record+1}=char(fread(fid,128,'char')');

end

fclose(fid);
 
% Find timestamps before break point
firstTSIdx=find(timeStamp<allTS(breakPt));
breakIdx=firstTSIdx(end);

timeStampFirst=timeStamp(1:breakIdx);
timeStampSecond=timeStamp(breakIdx+1:end);

numEventsRecFirst=length(timeStampFirst);
numEventsRecSecond=length(timeStampSecond);


% Now write the events files

% First file
fid=fopen(eventsFile,'w','l');
fwrite(fid,eventHeader,'char');

for record=0:numEventsRecFirst-1
  fwrite(fid,pktStart(record+1),'int16');
  fwrite(fid,pktId(record+1),'int16');
  fwrite(fid,pktDataSize(record+1),'int16');
  fwrite(fid,timeStampFirst(record+1),'int64');
  fwrite(fid,eventId(record+1),'int16');
  fwrite(fid,ttlValues(record+1),'int16');
  fwrite(fid,crc(record+1),'int16');
  fwrite(fid,dummy(record+1),'int32');
  fwrite(fid,extra(record*8+1:(record+1)*8),'int32');
  fwrite(fid,eventString{record+1}','char');

end

fclose(fid);


% Second file
fid=fopen(eventsFileSecond,'w','l');
fwrite(fid,eventHeader,'char');

for record=0:numEventsRecSecond-1
  fwrite(fid,pktStart(record+breakIdx+1),'int16');
  fwrite(fid,pktId(record+breakIdx+1),'int16');
  fwrite(fid,pktDataSize(record+breakIdx+1),'int16');
  fwrite(fid,timeStampSecond(record+1),'int64');
  fwrite(fid,eventId(record+breakIdx+1),'int16');
  fwrite(fid,ttlValues(record+breakIdx+1),'int16');
  fwrite(fid,crc(record+breakIdx+1),'int16');
  fwrite(fid,dummy(record+breakIdx+1),'int32');
  fwrite(fid,extra(record*8+breakIdx+1:(record+1)*8+breakIdx),'int32');
  fwrite(fid,eventString{record+breakIdx+1}','char');

end

fclose(fid);

