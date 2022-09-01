function rerefnlx(subjID,sessionID)
% rerefnlx - Rereferences neuralynx electrode signals based on macro
% micro arrangements
%
% NOTE - For now, this ignores the data collected on ref leads (RefA1 etc)
%
% Function:
% rerefnlx(subjID,sessionID)
%
% Input Args:
% 
%   subjID = 'UP011';                   % Subject ID
%   sessionID = 'UP011_17Jul07_1901';   % Session ID defined by date/time
%


% Create directories and filenames
elecDir = fullfile('/data/eeg', subjID, 'docs');
jackFile = fullfile('/data/eeg',subjID,'eeg.noreref.neuralynx',[sessionID '.jacksheet.txt']);
inDir = fullfile('/data/eeg',subjID,'eeg.noreref.neuralynx');
outDir = fullfile('/data/eeg',subjID,'eeg.reref.neuralynx');

% Create outDir
if ~exist(outDir,'dir')
  mkdir(outDir);
end

% Load in macro/micro grid data from electrode file
% Load in electrodes which is cell array of macro/micro grids
chdir(elecDir);
electrodes=[];
eval('electrodesNlx');

gridTypes=length(electrodes);

% Load in jacksheet data
[jackElec,csc]=textread(jackFile,'%d%s\n');

% Find CSC numbers
cscNum=[];
for i=1:length(csc)
     if csc{i}(1:3)=='CSC'
         cscNum=[cscNum;str2num(csc{i}(4:end))];
     end
end

% Now make new cell arrays, grid and elec, corresponding to new channel numbers

grid=cell(gridTypes);
elec=cell(gridTypes);


% First find where these cscs are in the macro and micro grids
for k=1:gridTypes
   for i=1:size(electrodes{k},1)
     idx=find(cscNum>=electrodes{k}(i,1) & cscNum<=electrodes{k}(i,2));
     if length(idx)>0
        grid{k}=[grid{k};jackElec(idx(1)) jackElec(idx(end))];
        elec{k}=[elec{k};jackElec(idx)];
     end
   end
end

weights=cell(gridTypes);

% Now get the weights for every grid type
for k=1:gridTypes

  weights{k}=ones(size(elec{k}));

  for i = 1:size(grid{k},1)
    % find the index of those leads which are in this grid
    idx = find(elec{k}>=grid{k}(i,1) & elec{k}<=grid{k}(i,2));
  
    if length(idx) > 0
      % weight them appropriately
      weights{k}(idx) = weights{k}(idx)/length(idx);
    end
  end

  % Normalize weights
  weights{k}=weights{k}/sum(weights{k});

end


% Create fileroot
fileroot=fullfile(inDir,sessionID);
fprintf('Processing %s...\n',fileroot);

% make a fake event to load data fromg gete
event = struct('eegfile',fileroot);

% get data info
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event);


% Now process data for different grid types

for k=1:gridTypes

  % Load leads for this gridType and calc avg
  fprintf('Calculating reference(%d): ',length(elec{k}));
  avg = [];
  for c = 1:length(elec{k})
    fprintf('%d ',c);
    teeg = gete(elec{k}(c),event,0);
    if isempty(avg)
      avg = (teeg{1}*weights{k}(c));
    else
      % Hack to fix the case when teeg{1} for last electrode contains 
      % smaller number of samples than avg
      if length(avg)>length(teeg{1})
        teeg{1}=[teeg{1},zeros(1,length(avg)-length(teeg{1}))];
      end
      avg = avg + (teeg{1}*weights{k}(c));
    end
  end
  fprintf('\n');
  
  % Load all leads for this gridType, apply avg, and save to new file
  fprintf('Saving rereferenced channels(%d): ',length(elec{k}));
  for c = 1:length(elec{k})
    fprintf('%d ',c);
    % load it
    teeg = gete(elec{k}(c),event,0);
    
    % Hack to fix the case when teeg{1} for last electrode contains 
    % smaller number of samples than avg
    if length(avg)>length(teeg{1})
      teeg{1}=[teeg{1},zeros(1,length(avg)-length(teeg{1}))];
    end

    % apply avg and reverse gain
    teeg{1} = (teeg{1}-avg)./gain;
    
    % save it
    [fdir,fname] = fileparts(fileroot);
    filestem = fullfile(outDir,fname);
    chanfile = sprintf('%s.%03i', filestem, elec{k}(c));
    % open and write the file
    fid = fopen(chanfile,'wb','l');
    fwrite(fid,teeg{1},dataformat);
    fclose(fid);
  end
  fprintf('\n');
  
end

% copy the params if there
pfile = fullfile(fileparts(fileroot),'params.txt');
if exist(pfile,'file')
  % copy to new location
  copyfile(pfile,outDir);
end


% copy jacksheet file
copyfile(jackFile,outDir);
