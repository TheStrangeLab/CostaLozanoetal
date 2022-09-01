function [upstrokes,downstrokes]=load_nev(filename)

% [uptimes, downtimes]=load_nev(filename)
% This function reads in the specified Neuralynx .Nev file, and outputs the
% upstrokes and downstrokes from it.
%
% From our setup at UCLA, we found that the TTL value of -2 corresponds to an
% upstroke, and -4 corresponds to a downstroke.  This function makes that assumption!

%revised 8/12: sometimes upstrokes seem to be -1, and downs seem to be
%-3. changed code accrdingly

% revised 7/23 for Penn neuralynx system - now upstrokes seem to be 1 and
% downstrokes seem to be 0.  there are occasional ttlvalues of 3, but we'll
% ignore these as they don't seem to be related to the actual sync pulses

% revised 3/24/2008 for UCLA sometimes upstrokes seem to be -13 and downstrokes
% -5. Also added new output for unclassified ttlValues.

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
  if ttlValue== -2 || ttlValue==-1 || ttlValue==1 || ttlValue == -13
    upstrokes=[upstrokes timeStamp];
  elseif ttlValue==-4 || ttlValue==-3 || ttlValue==0 || ttlValue == -5
    downstrokes=[downstrokes timeStamp];
  else
    disp(sprintf('%d is unclassified ttlValue!',ttlValue));
  end
  
end

fclose(fid);

