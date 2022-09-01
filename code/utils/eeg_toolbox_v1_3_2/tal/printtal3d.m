function printtal3d(res,filename)
%PRINTTAL3D - Prepare and print brain from tal3d res.
%
% FUNCTION:
%   printtal3d(res,filename);
%
% INPUT ARGS:
%   res % res returned from call to tal3d
%   filename = 'brain.png'; % filename to print
%
%

if ~exist('filename','var')
  filename = [];
end

if 0
  load RAW_coords.txt
  cxyz = RAW_coords;
  clear RAW_coords;
  %v3d = [-90,0];
end

% pick an alphacolor
acolor = [1 1 0];
%acolorint = uint8(fix(acolor.*double(intmax('uint8')))+1);
acolorint = [255 255 0];

axis off

% set the background to yellow for setting alpha later
set(gcf,'color',acolor,'InvertHardcopy','off')

% remove electrodes
set(res.hElec,'visible','off');

% print to temp file
fprintf('Saving brain...\n')
brainfile = [tempname '.png'];
print('-dpng','-r300',brainfile);

% replace electrodes and remove brain
set(res.hElec,'visible','on');
set(res.hBrain,'facealpha',0);

% print to temp file
fprintf('Saving electrodes...\n');
elecfile = [tempname '.png'];
print('-dpng','-r300',elecfile);

% load in the files
im_brain = imread(brainfile);
im_elecs = imread(elecfile);
sz = size(im_brain);

% set the original image back to normal
set(gcf,'color',[1 1 1],'InvertHardcopy','on')
set(res.hBrain,'facealpha',1);

% load the images onto surfaces
fprintf('Plotting brain and electrodes...\n');

figure
%clf
% set background back to white
%set(gcf,'color',[1 1 1],'InvertHardcopy','on')

% plot the brain
h1 = imagesc(im_brain);

% find where the image is the alphacolor and set alphadata
ind = find(~(im_brain(:,:,1)==acolorint(1) & im_brain(:,:,2)==acolorint(2) & im_brain(:,:,3)==acolorint(3))); 
adat = zeros(sz(1),sz(2));
adat(ind) = 1;
set(h1,'alphadata',adat)

% plot the elecs
hold on
h2 = imagesc(im_elecs);

% find where the image is the alphacolor and set alphadata
ind = find(~(im_elecs(:,:,1)==acolorint(1) & im_elecs(:,:,2)==acolorint(2) & im_elecs(:,:,3)==acolorint(3))); 
adat = zeros(sz(1),sz(2));
adat(ind) = 1;
set(h2,'alphadata',adat)

axis tight
axis off

hold off


% delete the temp files
delete(brainfile);
delete(elecfile);

% see if print to file
if ~isempty(filename)
  fprintf('Printing to file...');
  print('-dpng','-r300',filename);
end

fprintf('Done!\n');



