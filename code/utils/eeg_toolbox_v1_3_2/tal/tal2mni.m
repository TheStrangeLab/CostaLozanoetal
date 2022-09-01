function outpoints = tal2mni(inpoints)

% TAL2MNI - Talairach to MNI coordinates
% 
% outpoints = tal2mni(inpoints)
% 
% inpoints  - Nx3 or 3xN matrix of coordinates
%             (N being the number of points)
% 
% outpoints - the coordinate matrix with MNI points
% 
% See also, MNI2TAL & the best guess discussion at
% http://www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html
% 

% $Revision: 1.2 $ $Date: 2005/11/01 17:47:07 $

% Licence:  GNU GPL, no express or implied warranties
% Matthew Brett 2/2/01, matthew.brett@mrc-cbu.cam.ac.uk
% modified 02/2003, Darren.Weber_at_radiology.ucsf.edu
%                   - swapped inv() for slash equivalent
%                   - removed dependence on spm_matrix
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimdim = find(size(inpoints) == 3);
if isempty(dimdim),
    error('input must be a Nx3 or 3xN matrix')
end
if dimdim == 2 | length(dimdim) == 2,
    inpoints = inpoints';
end

% Transformation matrices, different zooms above/below AC
M2T = mni2tal_matrix;

inpoints = [inpoints; ones(1, size(inpoints, 2))];

tmp = inpoints(3,:) < 0;  % 1 if below AC

inpoints(:,  tmp) = (M2T.rotn * M2T.downZ) \ inpoints(:,  tmp);
inpoints(:, ~tmp) = (M2T.rotn * M2T.upZ  ) \ inpoints(:, ~tmp);

outpoints = inpoints(1:3, :);
if dimdim == 2 | length(dimdim) == 2,
    outpoints = outpoints';
end

return
