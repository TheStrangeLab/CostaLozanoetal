function res = region3d(viewAZEL,regioncolors,defaultcolor)
%REGION3D - Topographical map of colored 3D brain surface.
%
% REGIONCOLORS is an Nx2 cell matrix of location lookup parameters and
% the color assigned to them (see the example row below).
%
% If you leave REGIONCOLORS blank and set DEFAULTCOLOR to a single
% value between 1 and 5 corresponding to the database indexes, it
% will plot random colors for each lookup.  
%
% To see the five possible location lookups, load the file
% mri_cortical_surface.mat and examine the structure loc_lookup.
% Use the following names to refer to the lookups in your
% REGIONCOLORS specification:
%
% 1: hemisphere
% 2: lobe- Frontal, temporal, etc ...
% 3: label- Essentially the gyri
% 4: type- Gray matter, white matter, CSF
% 5: brodmann- Brodmann areas
% 
% FUNCTION:
%   res = region3d(viewAZEL,regioncolors,defaultcolor)
%
% INPUT ARGS:
%   viewAZEL = [90 0]; % Will be a right lateral view
%   regioncolors = {'strfound(hemisphere,''Right'')&strcmp(label,''Middle Temporal Gyrus'')',[1 0 0]}; 
%   defaultcolor = [.5 .5 .5]; % Default brain color
%
%
% OUTPUT ARGS:
%   res.hLight - handle to the camera light.
%   res.hBrain - handle to the brain.
%
%

conversion = {'hemisphere','loc_lookup{1}';
              'lobe','loc_lookup{2}';
              'label','loc_lookup{3}';
              'type','loc_lookup{4}';
              'brodmann','loc_lookup{5}';};


% process input vars
if ~exist('viewAZEL','var') | isempty(viewAZEL)
  viewAZEL = [90,0];
end
% process input vars
if ~exist('regioncolors','var')
  regioncolors = [];
end
% process input vars
if ~exist('defaultcolor','var') | isempty(defaultcolor)
  defaultcolor = [.5 .5 .5];
end


% set up return structure
res.hLight = [];
res.hBrain = [];

% read in the faces and vertices of the surface
% get path to pictures
picpath = fileparts(which('region3d'));
picfile = fullfile(picpath,'mni_cortical_surface.mat');
load(picfile);

% ignore long distance locs
vloc(vdist>10) = 0;

% calc the cdata
if length(defaultcolor) == 1
  % put random data in the specified lookup column
  l = defaultcolor;
  cdata = zeros(size(v,1),3);
  for g = 1:length(loc_lookup{l})
    ind = find(vloc(:,l)==g);
    cdata(ind,:) = ones(length(ind),1)*rand(1,3);
  end  
else
  % use the regioncolors
  cdata = ones(size(v,1),1)*defaultcolor;
  if ~isempty(regioncolors)
    %fprintf('Processing regions...');
    for i = 1:size(regioncolors,1)
      ind = findloc(regioncolors{i,1},loc_lookup,vloc);
      if ~isempty(ind)
        cdata(ind,:) = ones(length(ind),1)*regioncolors{i,2};
      end
    end
  end
end

%fprintf('Done\n');

%l = 3;
%g =5;
%loc_lookup{l}{g}
%cdata = ones(size(v,1),3)*.5;
%ind = find(vloc(:,l)==g);
%cdata(ind,:) = ones(length(ind),1)*[1 0 0];


% draw the patches
hs = patch('faces',f,'vertices',v,'EdgeColor','none','FaceColor','interp','FaceVertexCData',cdata);
daspect([1 1 1]);
view(viewAZEL)
res.hLight = camlight;
set(res.hLight,'Color',[.8 .8 .8],'Style','infinite');
lighting phong
axis off

% save the handle
res.hBrain = hs;





%%%%%%%%%%%%%%%%%
function ind = findloc(filterstr,loc_lookup,vloc)
%
%

conversion = {'hemisphere','loc_lookup{1}';
              'lobe','loc_lookup{2}';
              'label','loc_lookup{3}';
              'type','loc_lookup{4}';
              'brodmann','loc_lookup{5}';};
              

% loop over conversions and replace
for i = 1:size(conversion,1)
  filterstr = strrep(filterstr,conversion{i,1},conversion{i,2});
end

% find the & and |
ands = strfind(filterstr,'&');
ors = strfind(filterstr,'|');
both = sort([ands ors]);

% eval and replace str comp with number comparisons
newfilterstr = '';
numexp = length(both)+1;
istart = 1;
for i = 1:numexp
  if i==numexp
    iend = length(filterstr);
  else
    iend = both(i)-1;
  end
  
  % get the compare indexes
  evalstr = filterstr(istart:iend);
  cind = eval(evalstr);
  cind = find(cind);

  if isempty(cind)
    % no areas found, so no color
    newfilterstr = '0';
    break;
  end
  
  % replace the string with the correct comparison
  newstr = '[';
  column = '';
  for c = 1:size(conversion,1)
    if strfound(evalstr,conversion{c,2})
      % now we have the correct column
      column = ['vloc(:,' num2str(c) ')'];
    end
  end
  for c=1:length(cind)
    newstr = [newstr column '==' num2str(cind(c))];
    
    if c ~= length(cind)
      % not at end, so must add |
      newstr = [newstr ' | '];
    end
  end
  
  % close the newstr
  newstr = [newstr ']'];
  
  % replace the chunk in filterstr
  newfilterstr = [newfilterstr newstr ];
  
  % see if must add conjunction
  if i~=numexp
    newfilterstr = [newfilterstr filterstr(both(i))];
  end
  
  % set the new start
  istart = iend+2;  
end

% eval that action
vals = eval(newfilterstr);
ind = find(vals);

