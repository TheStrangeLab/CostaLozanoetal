function mat = getmatrix(filename,col)
%GETMATRIX - Loads a n-dim matrix of singles from a file.
%
% Loads a ndim dimensional matrix from a binary file saved with
% putmatrix.  The PUTMATRIX / GETMATRIX combo of functions is a
% good way to save multidimensional matrixes to files without
% taking up loads of space.
%
% If you are dealing with very large matrixes, then you can specify
% a single column of data from the file to load with the COL
% parameter.  With multidimensional matrixes, this specifies the
% final dimension and will return all other dimensions.
%
% FUNCTION:
%   mat = getmatrix(filename,col)
%
% INPUT ARGS:
%   filename = 'recalled/energy/pow_01.dat';
%   col = 10; % (Optional) Column of data to read in
%
% OUPUT ARGS:
%   mat - the matrix of data returned
%

if nargin < 2
  col = 0;
end

% set the single size
single_size = 4;

% open the file
fid = fopen(filename,'rbl');

% get the number of dimensions and then dimensions
ndim = fread(fid,1,'single');
dims = zeros(1,ndim);
for i = 1:ndim
  dims(i) = fread(fid,1,'single');
end


% see if we are skiping to a certain column
if col > 0 & ndim > 0
  stat = fseek(fid,(col-1)*single_size*prod(dims(1:ndim-1)),'cof');
  if stat < 0
    % error
    fprintf('%s\n',ferror(fid));
    fclose(fid);
    return
  end
  
  % read it in
  mat = fread(fid,prod(dims(1:ndim-1)),'single');
  
  % reshape it
  mat = single(reshape(mat,dims(1:ndim-1)));
else

  % read it in, yo
  mat = fread(fid,inf,'single');
  
  mat = single(reshape(mat,dims));
  
end


% close it up
fclose(fid);

