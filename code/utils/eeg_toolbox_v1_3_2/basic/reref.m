function reref(fileroots,grids,outdir,taldir)
%REREF - Rereference EEG recording with weights based on grids.
%
% Requires: leads.txt and good_leads.txt
%
% FUNCTION:
%   reref(fileroots,grids,outdir,taldir)
%
% INPUT ARGS:
%   fileroots = events;     % Events struct or cell array of EEG data files.
%   grids = [1,8;           % Each row is beginning and end of each grid.
%            9,24;
%            25,40;
%            41,48;
%            49,56;
%            57,64;];
%   outdir = './eeg.reref'; % Dir where new files will be saved.
%   taldir = './tal' or {all,good};       % Where to find leads.txt and good_leads.txt
%
%
%


if nargin < 4
  taldir = './tal';
  if nargin < 3
    outdir = './eeg.reref';
  end
end

% create dir
if ~exist(outdir,'dir')
  mkdir(outdir);
end

% load the leads
if isstr(taldir)
  % load from file
  alfile = fullfile(taldir,'leads.txt');
  glfile = fullfile(taldir,'good_leads.txt');
  al = getleads(alfile);
  gl = getleads(glfile);
else
  % get from cell
  al = taldir{1};
  gl = taldir{2};
end

% get the weights
weights = ones(size(gl));
for i = 1:size(grids,1)
  % find the index of those good leads which are in this grid
  idx = find(gl>=grids(i,1) & gl<=grids(i,2));
  
  if length(idx) > 0
    % weight them appropriately
    weights(idx) = weights(idx)/length(idx);
  end
end

% normalize up front
weights=weights/sum(weights); 

% process the fileroots
if isstruct(fileroots)
  % is events struct, so pull out unique file roots
  events = fileroots;
  fileroots = unique(getStructField(events,'eegfile','~strcmp(eegfile,'''')'));
end

% loop over the file roots
for f = 1:length(fileroots)
  fprintf('Processing %s...\n',fileroots{f});
  % make a fake event to load data fromg gete
  event = struct('eegfile',fileroots{f});
  
  % get data info
  [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event);

  % Load good leads and calc avg
  fprintf('Calculating reference(%d): ',length(gl));
  avg = [];
  for c = 1:length(gl)
    fprintf('%d ',c);
    teeg = gete(gl(c),event,0);
    if isempty(avg)
      avg = (teeg{1}*weights(c));
    else
      avg = avg + (teeg{1}*weights(c));
    end
  end
  fprintf('\n');
  
  % Load all leads, apply avg, and save to new file
  fprintf('Saving rereferenced channels(%d): ',length(al));
  for c = 1:length(al)
    fprintf('%d ',c);
    % load it
    teeg = gete(al(c),event,0);
    
    % apply avg and reverse gain
    teeg{1} = (teeg{1}-avg)./gain;
    
    % save it
    [fdir,fname] = fileparts(fileroots{f});
    filestem = fullfile(outdir,fname);
    chanfile = sprintf('%s.%03i', filestem, al(c));
    % open and write the file
    fid = fopen(chanfile,'wb','l');
    fwrite(fid,teeg{1},dataformat);
    fclose(fid);
  end
  fprintf('\n');
  
  % copy the params if there
  pfile = fullfile(fileparts(fileroots{f}),'params.txt');
  if exist(pfile,'file')
    % copy to new location
    copyfile(pfile,outdir);
  end
  
end





