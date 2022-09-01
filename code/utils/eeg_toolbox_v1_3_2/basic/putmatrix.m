function putmatrix(OutFileName,OutData)
%PUTMATRIX - Saves an n-dim matrix of singles to a file.
%
% Saves an ndim dimensional matrix from a binary file.  The 
% PUTMATRIX / GETMATRIX combo of functions is a good way to 
% save multidimensional matrixes to files without taking up 
% loads of space.
%
% FUNCTION:
%   putmatrix(OutFileName,OutData)
%
% INPUT ARGS:
%   OutFileName = 'pow_01.dat';  % Filename to save
%   OutData = recpow; % Data to save to file
%

% save to file
fid = fopen(OutFileName,'wbl');

% write out the number of dimensions
fwrite(fid,length(size(OutData)),'single');

% write out the value of the dimensions
fwrite(fid,size(OutData),'single');

% write out the data
fwrite(fid,OutData,'single');

fclose(fid);
