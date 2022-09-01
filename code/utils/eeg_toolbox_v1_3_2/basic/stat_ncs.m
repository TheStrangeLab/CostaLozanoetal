function info=stat_ncs(filename)
%function info=stat_ncs(filename)
%
%This function returns stats about a .Ncs file in a struct
% info.numRecords: number of records (512 samples per record).
% info.numSamples: number of samples in file (eqiv. to numRecords*512)
% info.sampleTSDiff: time in microseconds between samples
% info.sampleRate: sampling rate
% info.firstSampleTime: cheetah time of FIRST sample
% info.lastSampleTime: cheetah time of LAST sample
% info.ADGain: A/D gain
% info.AMPGain: AMP gain
% info.bitVolts: volts / bit


bytesPerRecord=1024+20;
headerBytes=16*1024;
samplesPerRecord=512;

d=dir(filename);
info.numRecords=(d.bytes-headerBytes)/bytesPerRecord;

if rem(info.numRecords,1)~=0
  error('length of file not an even number of blocks');
end

fid=fopen(filename,'r','l');

%read the header into a string
headerString=fread(fid,headerBytes,'*char')';

%then read the first two timestamps
%ts=fread(fid,2,'int64',bytesPerRecord-8); %8 is the size of an
%int64
ts=fread(fid,inf,'int64',bytesPerRecord-8); %8 is the size of an int64
fseek(fid,headerBytes+16,'bof');
numValid=fread(fid,'int32',bytesPerRecord-4); %8 is the size of an int32

if any(numValid~=samplesPerRecord)
  error('some non-full blocks!');
end

% Recordtsdiff=diff(ts);
% sampleTSDiff=recordTSDiff/samplesPerRecord;
% info.sampleRate=1/(1e-6*sampleTSDiff); %cheetah times are in microseconds

info.headerSampleRate=regexp_grab_token(headerString,'SamplingFrequency\s+([\d\\.]+)');

allDiffs=diff(ts);
goodDiffs=allDiffs(allDiffs<2*median(allDiffs));
sampleTSDiff=mean(goodDiffs)/samplesPerRecord;
info.actualSampleRate=1e6/sampleTSDiff;

info.firstSampleTime=ts(1);
lastBlockStart=ts(end);

%find the last time stamp, by getting the time of the last RECORD, then
%computing what the last time in that BLOCK should be.

info.lastSampleTime=lastBlockStart+(samplesPerRecord-1)*sampleTSDiff;

fclose(fid);


info.bitVolts=regexp_grab_token(headerString,'ADBitVolts\s+([\d\\.]+)');
%info.ADgain=regexp_grab_token(headerString,'ADGain\s+([\d\\.]+)');
%info.AMPgain=regexp_grab_token(headerString,'AmpGain\s+([\d\\.]+)');

info.numPotentialSamples=round((info.lastSampleTime-info.firstSampleTime)/sampleTSDiff)+1;

info.numActualSamples=info.numRecords*samplesPerRecord;

dateString=regexp_grab_token(headerString,'Opened \(m/d/y\): ([\d/]+)');
pat='Opened \(m/d/y\): [\d/]+\s+At Time: ([\d:\.])+';
timeString=regexp_grab_token(headerString,pat);

info.fileStartTime=[dateString ' ' timeString];

info.sampleTSDiff=sampleTSDiff;


function n=regexp_grab_token(string,pat)
%takes a regexp that captures one token, and returns it as a number (if it
%can). 
%find_num_in_str(string,pat)
[s,f,t]=regexp(string,pat);
if length(t)==0
  n=[];
  return;
end

token=string(t{1}(1):t{1}(2));

%n=str2num(token);
n=str2double(token);

if isnan(n) %if we couldn't return it as a number
  n=token; %then return a string. better than nothing, i guess.
end

