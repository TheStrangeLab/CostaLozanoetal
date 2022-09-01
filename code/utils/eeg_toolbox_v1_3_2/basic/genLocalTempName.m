function tmpFile = genLocalTempName(filename)
%GENLOCALTEMPNAME - Generate a temp file based on filename.
%
% This function will generate a filename for use as a temporary
% file that is local to the machine you are running on.
%
% FUNCTION:
%   tmpFile = genLocalTempName(filename)
%
%
%



% get the hostname
[s,hostname] = system('hostname');
hostname = strtrim(hostname);

% get a tempfile
[path,tfile] = fileparts(tempname);

% create the tempfile
tmpFile = [filename '.' hostname '.' tfile '.temp'];


