function convertFile_mni2tal(coordsfile)
%CONVERTFILE_MNI2TAL - Convert RAW_coords file from mni to tal.
%
% The new Talairach-coordinates file will be named RAW_coords.txt. The
% backup MNI-coordinates file will be renamed with a .mni extension if
% it's not already.
%
% FUNCTION:
%   convertFile_mni2tal(coordsfile)
%
% INPUT ARGS:
%   coordsfile = 'RAW_coords.txt.mni';  % File to convert
%
%
%



% read in the file
[c,x,y,z] = textread(coordsfile,'%d%n%n%n');

mni_coords = [x y z];

% convert the coords
tal_coords = mni2tal(mni_coords);

% backup the old file
if strcmp(coordsfile(end-3:end),'.mni')
  backupfile = coordsfile;
  fprintf('MNI backup file %s exists\n',backupfile);
  coordsfile = coordsfile(1:end-4);
else
  backupfile = [coordsfile '.mni'];
  fprintf('Backing up MNI file to %s\n',backupfile);
  movefile(coordsfile,backupfile);
end

% write out new file
outdata = [c tal_coords];
fid = fopen(coordsfile,'w');
fprintf(fid,'%d\t%g\t%g\t%g\n',outdata');
fclose(fid);

fprintf('Wrote new Tal file to %s\n',coordsfile);
