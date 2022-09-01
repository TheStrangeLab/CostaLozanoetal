function view3d(res,v)
%VEIW3D - Rotate a 3d image to new view, preserving lighting.
%
% This function takes in the struct of handles returned from a call
% to either tal3d or region3d and will rotate the image to the
% specified view, restoring the light parameters to be from the
% top-right and infinte style.
%
% FUNCTION:
%   view3d(res,v);
%
% INPUT ARGS:
%   res- struct of handles returned from call to region3d or tal3d
%   v = [90,0];  % view to set to (see view function for info).
%                % some good values are:
%                %  inferior = [-180 -90];
%                %  right sag. = [90 0];
%                %  left sag. = [-90 0];
%

% set the view
view(v);

if isfield(res,'hLight')
  % update the light
  camlight(res.hLight,'RIGHT');
  set(res.hLight,'style','infinite');
end



