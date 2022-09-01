function var_to_save = loadVar(filename)
%LOADEVENTS - Load a variable from a file.
%
% Use this function to load a variable from a file and set
% it to a specified variable.  The .mat file must have been saved
% using saveVar.
%
% FUNCTION:
%   varout = loadVar(filename)
%
% INPUT ARGS:
%   filename = 'events/events.mat';  % .mat file an events struct
%
% OUTPUT ARGS:
%   varout - The var that was saved to the file.
%


% load the mat file
load(filename)

% make sure the var exists
if ~exist('var_to_save','var')
  var_to_save = [];
  warning('File was not saved with saveVar.');
end





