function created = touchfile(filename,pausetime)
%TOUCHFILE - Check if file exists, touching if not.
%
% This function helps to keep multiple procs from working on the
% same file when on a cluster.  It will create a temporary file and
% wait two seconds for another process to try and claim the file.
% If no other process claims the file, then it will touch the file
% and clean up.
%
% FUNCTION:
%   created = touchfile(filename,pausetime)
%
% INPUT ARGS:
%   filename- The file to test/touch.
%   pausetime- Optional pause time in seconds. (Default is 2).
%
% OUTPUT ARGS:
%   created- 1 if we created the file and are responsible for it
%
%

if ~exist('pausetime','var')
  pausetime = 2;
end

% test name
checkname = [filename '.check'];

% first see if file exists
if exist(filename,'file') | exist(checkname,'file')
  % some proc is already wanting to start
  created = 0;
  return;
end

% touch the check file
system(['touch ' checkname ' ; sync']);

% pause to allow for other procs to have a chance
pause((pausetime + rand*2 - 1));

% see if the file still does not exist
if exist(filename,'file')
  % some proc is already wanting to start
  if exist(checkname,'file')
    system(['rm -f ' checkname '; sync']);
  end
  created = 0;
  return;
end

% touch the file
system(['touch ' filename ' ; sync']);
created = 1;

% clean up if necessary
if exist(checkname,'file')
  system(['rm -f ' checkname '; sync']);
end

