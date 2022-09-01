function dat = look(filestem, e1, e2)  
%LOOK - Load and look at EEG data.
%
% Use LOOK to load all the data from an eeg file and plot it.  It
% will return the data for further analysis.  You can use it in
% three ways:
%
% 1) >> look('CH012_dat.001');  % load from the file
% 2) >> look('CH012_dat',1);    % appends the channel to the root
% 3) >> look('CH012_dat',1,2);  % loads two channels and takes diff
%
% FUNCTION:
%   dat = look( filestem, [e1 [, e2]] )
%
% INPUT ARGS:
%   filestem = 'CH012_dat.001';  % filestem or file to load
%   e1 = 2;   % (Optional) Lead to append to stem
%   e2 = 3;   % (Optional) Bipolar lead to subtrace from e1
%
% OUTPUT ARGS:
%   dat - Data samples from the file
%

% 3/6/04 PBS - Looks for params.txt file to load rate, format, and gain.
% 3/23/04 JJ - Put in code to handle the filestem containing a trailing '.'


if filestem(end)=='.'
  filestem=filestem(1:end-1);
end

% get the data format
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(fileparts(filestem));

% read in the data
switch nargin
 case 1,
  %f = openfile( filestem, 'r', 'l');
  f = fopen( filestem, 'r','l');
  dat =  fread(f, inf, dataformat);
 case 2,
  %f = openfile( sprintf('%s.%.3i', filestem, e1), 'r', 'l');
  f = fopen( sprintf('%s.%.3i', filestem, e1), 'r','l');
  dat =  fread(f, inf, dataformat);
 case 3,
  %f = openfile( sprintf('%s.%.3i', filestem, e1), 'r', 'l');
  f = fopen( sprintf('%s.%.3i', filestem, e1), 'r','l');
  dat =  fread(f, inf, dataformat);
  fclose(f);
  %f = openfile( sprintf('%s.%.3i', filestem, e2), 'r', 'l');
  f = fopen( sprintf('%s.%.3i', filestem, e2), 'r','l');
  dat =  [dat fread(f, inf, dataformat)];
  dat = dat(:,2) - dat(:,1);
 otherwise,
  help look;
  return;
end

% apply the gain
dat = dat.*gain;

% plot it
plot(dat);
if nargin == 1
  title( filestem );
else
  title( sprintf('%s %i', filestem, e1) );
end
fclose(f);


 
if nargout == 0
  clear dat;
end
%figure(fh);
