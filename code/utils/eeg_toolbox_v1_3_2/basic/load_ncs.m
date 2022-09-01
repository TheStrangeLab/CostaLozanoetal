function [data,info,allTS]=load_ncs(filename)
%[data,info]=load_ncs(filename)
%
%this function returns the raw 2 byte samples inside this .Ncs file
%this value should be multipled by the info.ADGain (from stat_ncs) to get the
%actual voltage value.

samplesPerRecord=512;
bytesPerRecord=1024+20;
headerBytes=16*1024;

info=stat_ncs(filename);
disp(sprintf('difference between samples is %g microseconds',info.sampleTSDiff));
disp(sprintf('header sampling rate is %g',info.headerSampleRate));
disp(sprintf('actual sampling rate is %g',info.actualSampleRate));

if info.numPotentialSamples~=info.numActualSamples
  warning('missing some samples...');
end

data=zeros(1,info.numPotentialSamples,'int16'); %preallocation
%data=zeros(1,info.numActualSamples,'int16'); %preallocation
fid=fopen(filename,'r','l');

fseek(fid,headerBytes,'bof'); 

prevTS=nan;
allTS=zeros(1,info.numRecords);

for record=0:info.numRecords-1
  ts=fread(fid,1,'int64');

  tmp=fread(fid,3,'int32'); %read in 3 in one shot for speed
  numValidSamp=tmp(3);
  
  
  thisRecordData=fread(fid,512,'*int16');
  firstIndex=round((ts-info.firstSampleTime)/info.sampleTSDiff);
  i=1:numValidSamp;
  data(firstIndex+i)=thisRecordData(i);

  
  if ~isnan(prevTS)& abs(ts-(prevTS+info.sampleTSDiff*512))>10 
    %if difference is greater than 10 microseconds
    %note: 10 microseconds is totally arbitrary... someone should
    %think about this more.  i guess this shit is due to rounding,
    %as indicated by some .pdf in the neuralynx documentation, but
    %this does seem like a rather large difference....
    
    warning('timestamp difference has changed... yuck');
  elseif (numValidSamp~=512)
    error('512 samples WERE NOT read. WTF');
  end
  
  prevTS=ts;
  allTS(record+1)=ts;
end

fclose(fid);
