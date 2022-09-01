function p=eegparams(field,paramdir)
%EEGPARAMS - Get a subject specific eeg parameter from the params.txt file.
% 
% If paramdir is not specified, the function looks in the 'docs/'
% directory for the params.txt file.
%
% The params.txt file can contain many types of parameters and will
% evaluate them as one per line.  These are examples:
%
% Channels 1:64
% samplerate 256
% subj 'BR015'
%
% FUNCTION:
%   p = eegparams(field,paramdir)
%
% INPUT ARGS:
%   field = 'samplerate';        % Field to retrieve
%   paramdir = '~/eeg/012/dat';  % Dir where to find params.txt
%
% OUTPUT ARGS:
%   p- the parameter, evaluated with eval()
%

% VERSION HISTORY:
%

if nargin < 2
  paramfile = 'docs/params.txt';
else
  paramfile = [paramdir '/params.txt'];
end


if(strcmp(field,'samplerate')) p=256.03; % default, BioLogic standard
elseif(strcmp(field,'gain')) p=1; % default, BioLogic standard
else p=[]; end; % returns an empty matrix if not a known field and not
                % found in the file

		
in=fopen(paramfile,'rt');
if(in~=-1)
  done=0;
  f=0;
  while( (~isempty(f)) & ~done)
    f=fscanf(in,'%s',1); 
    %v=fscanf(in,'%f',1);
    v = fgetl(in);
    if(strcmp(f,field)) 
      if isstr(v)
	v = eval(v);
      end
      
      p=v; 
      done=1; 
    end;

  end % while not done
  fclose(in);
end
% default is to return the default value assigned at the top



% this way is slooow

% % load the params
% if exist(paramfile,'file')
%   [fname,param] = textread(paramfile,'%s%[^\n]');
  
%   ind = find(strcmp(fname,field));
%   if ~isempty(ind)
%     % set the value
%     p = eval(param{ind(1)});
%   end
% end
