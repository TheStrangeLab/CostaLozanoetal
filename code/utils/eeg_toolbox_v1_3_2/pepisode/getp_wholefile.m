function uvec=getp_wholefile(dataDir,chan)
%Input:
%  dataDir- directory containing pepisode data (may need to start with
%  'pepisode' for compatibility with calcPepisode.m)
%  chan- channel number to look at
%
%Output:
%  uvec - boolean matrix of pepisode values for entire session (has dimensions
%  frequency X sample)


dataFile=fullfile(dataDir,sprintf('chan_union.%.03i.gz',chan));
uncompressedName=[tempname sprintf('%.f',1e3*rand)];

system(sprintf('gunzip -c %s > %s',dataFile,uncompressedName));


fid=fopen(uncompressedName,'r');

nrows = fread(fid,1,'int');
ncols = fread(fid,1,'int');
uvec=fread(fid,[nrows ncols],'char')==1;

fclose(fid);
delete(uncompressedName);
