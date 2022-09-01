function new_coords=get_tal_montage(montagename)
% new_coords=get_tal_montage(montagename)
%
% This function returns the electrode coordinates for the electrode
% numbers listed in the passed montage file. The montage file for the
% montagename 'inf' should be named tal/inf.montage; same for lsag
% and rsag.

mfname=sprintf('tal/%s.montage',montagename);
in=fopen(mfname,'r');
if(in==-1)
  fprintf(1,'No montage file %s.\n',mfname); new_coords=[];
else
  selectthese=fscanf(in,'%i',inf); fclose(in);
  if exist('./RAW_coords.txt','file') | exist('tal/RAW_coords.txt','file')
    new_coords=get_tal_coords('RAW_coords.txt',selectthese);
  else
    new_coords=get_tal_coords('raw_coords.txt',selectthese);
  end
end
