function dat = readbinary(filename,size,format)
%READBINARY - Wrapper function to read in data from a binary file.
%
%
% FUNCTION:
%   dat = readbinary(filename,size,format);
%
%



[fid,message] = fopen(filename,'rb');
if fid == -1
  disp(message)
  return
end

dat = fread(fid,size,format);

fclose(fid);

