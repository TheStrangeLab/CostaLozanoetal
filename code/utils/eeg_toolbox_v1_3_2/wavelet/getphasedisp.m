function [p_ray,p_ray_boot,disp_diff,disp_diff_boot,phase] = getphasedisp(chan,events,DurationMS,OffsetMS,BufferMS,ShiftMS,varargin)
%GETPHASEDISP - Test non-uniformity and phase dispersion for a set of events.
% 
% Calculates the wavelet phase as a funtion of time and frequency
% for the specified sets of events, calculates the rayleigh 
% uniformity statistic for each point for each set of events, then
% calculates the difference in dispersion between each pairwise set
% of events.
%
% Then you can optionally perform a bootstrap on all the tests by
% setting the ITERATIONS parameter.
%
% The Rayleigh bootstrap entails randomly shifting the signal of
% each event.  The dispersion ANOVA test randomly permutes the events
% from the two groups.
%
% FUNCTION:
%   [p_ray,p_ray_boot,disp_diff,disp_diff_boot,phase] = getphasedisp(chan,events,DurationMS,OffsetMS,BufferMS,ShiftMS,varargin)
%
% INPUT ARGS:
%   chan = 2;     % Channel on which to perform the tests
%   events = {rec_events,not_events}; % Cell array of sets of
%                                     %   events structures
%   DurationMS = 2500;  % Duration over which to calculate
%   OffsetMS = -500;    % Offset from start of event
%   BufferMS = 1000;    % Buffer used for calculation
%   ShiftMS = 1000;     % +/- Shift used in the Rayleigh bootstrap
%   
%   OPTIONAL PARAMS:
%     'freqs', [2,4,8]  % Frequencies to analyze, defaults to
%                       %   eeganalparams('freqs')
%     'width', 6        % Width of the Morlet wavelet, defaults to
%                       %   eeganalparams('width')
%     'filtfreq'
%     'filttype'
%     'filtorder'
%     'iterations', 100 % perform a bootstrap
%
% OUTPUT ARGS:
%   p_ray- Cell array of p values from the rayleigh stat, one
%            time-frequency matrix per set of events
%   disp_diff- Cell array (One cell per event) of the results from
%            the disp_anova test on each pair of data.  
%   phase- The actual phase values for each event, freq, and time
%   
%   
% NOTE: Requires PhasePACK
%

% 2/10/05 - PBS - Change disp_anova sig values from p*dir to be
%                 either p or 1-p depending on the direction.
% 1/20/05 - PBS - Fixed factorial bug.
%
%fixed a bug with the filtfreq function

% events must be a cell array
if ~iscell(events)
  error('events must be a cell array.');
end

% set the defaults
freqs = eeganalparams('freqs');
width = eeganalparams('width');
resampledrate = [];
filtfreq = []; %changed by josh
filttype = 'stop';
filtorder = 1;
iterations = 0;

do_rayleigh = 1;
do_disp_anova = 1;

% process the varargs
i = 1;
while length(varargin)>0 & i<=length(varargin)
  switch lower(varargin{i})
   case 'freqs'
    i = i+1;
    freqs = varargin{i};
    i = i+1;
   case 'width'
    i = i+1;
    width = varargin{i};
    i = i+1;
   case 'resampledrate'
    i = i+1;
    resampledrate = varargin{i};
    i = i+1;
   case 'filtfreq'
    i = i+1;
    filtfreq = varargin{i};
    i = i+1;
   case 'filttype'
    i = i+1;
    filttype = varargin{i};
    i = i+1;
   case 'filtorder'
    i = i+1;
    filtorder = varargin{i};
    i = i+1;
   case 'iterations'
    i = i+1;
    iterations = varargin{i};
    i = i+1;
    
   case 'no_rayleigh'
    do_rayleigh = 0;
    i = i + 1;
   case 'no_disp_anova'
    do_disp_anova = 0;
    i = i + 1;
    
   otherwise
    error(['Error processing vararg: ' num2str(i)]);
  end
end

% get some parameters
rate = eegparams('samplerate',fileparts(events{1}(1).eegfile));
samplerate = rate;
if ~isempty(resampledrate)
  resampledrate = round(resampledrate);
  samplerate = resampledrate;
  rate = resampledrate;
end

% convert the durations to samples
duration = round((DurationMS)*rate/1000);
shift = round((ShiftMS)*rate/1000);
offset = round((OffsetMS-BufferMS)*rate/1000);
buffer = round((BufferMS)*rate/1000);

% get the min n for the events
for e = 1:length(events)
  min_n(e) = length(events{e});
end
min_n = min(min_n);

% num comparisons
numcomp = (length(events)*(length(events)-1))/2;

% get phase for each set of events, adding in the amount needed for shift
fprintf('Getting phase ...\n');
phase = cell(1,length(events));
for i = 1:length(events)
  phase{i} = getphasepow(chan,events{i},DurationMS+(2*ShiftMS),OffsetMS-ShiftMS,BufferMS,'filtfreq',filtfreq,'filttype',filttype,'filtorder',filtorder,'width',width,'freqs',freqs,'resampledrate',resampledrate);
end
fprintf('\n');


% allocate some space
p_ray = [];
p_ray_boot = [];
disp_diff = [];
disp_diff_boot = [];

if do_rayleigh
  p_ray = cell(2,length(events));
  p_ray_boot = cell(1,length(events));
  for i=1:length(events)
    p_ray{1,i} = single(zeros(length(freqs),duration));
    p_ray{2,i} = single(zeros(length(freqs),duration));

    if iterations > 0
      p_ray_boot{1,i} = single(zeros(iterations,length(freqs),duration));
    end
  end
end

if do_disp_anova
  disp_diff = cell(2,numcomp);
  disp_diff_boot = cell(2,numcomp);
  for i=1:numcomp
    disp_diff{1,i} = single(zeros(length(freqs),duration));
    disp_diff{2,i}  = single(zeros(length(freqs),duration));
    
    if iterations > 0
      disp_diff_boot{1,i} = single(zeros(iterations,length(freqs),duration));
      disp_diff_boot{2,i}  = single(zeros(iterations,length(freqs),duration));
    end
  end
end


% loop and calc phase dispersion and diff over time
fprintf('Statistics...\nFreqs(%d): ',length(freqs));
for f = 1:length(freqs)
  fprintf('%d ',f);

  for t= (1+shift):(duration+shift)  % adjusted to not do the shift area

    if do_rayleigh
      % perform uniformity test for each
      for i = 1:length(events)
	[p_ray{1,i}(f,t-shift),p_ray{2,i}(f,t-shift)] = rayleigh_alt(double(phase{i}(:,f,t)),min_n);
      end
    end
    
    if do_disp_anova
      % perform disp_anova test on all pairwise combinations
      count = 0;
      for i = 1:length(events)-1
	for j = i+1:length(events)
	  count = count+1;
	  
	  [disp_diff{1,count}(f,t-shift),disp_diff{2,count}(f,t-shift),direction] ...
	      = disp_anova(double(phase{i}(:,f,t)),double(phase{j}(:,f,t)));
	  
	  %disp_diff{1,count}(f,t-shift) = single(double(disp_diff{1,count}(f,t-shift))*direction);
          disp_diff{1,count}(f,t-shift) = single( (1 + double(disp_diff{1,count}(f,t-shift))*direction) + ((-direction - 1)/2) );
	end
      end
    end
  end
  
  % see if do a bootstrap
  if iterations > 0
    fprintf('\nIterations(%d):\n',iterations);
    
    for loop=1:iterations
      fprintf('%d ',loop);

      if do_rayleigh
	% loop over groups of events
	for i = 1:length(events)
	  % randomize time shift for each event
	  tstart = round((2*shift*rand(1,size(phase{i},1)))-.5) + 1;
	  
	  ptemp = squeeze(phase{i}(:,f,:));
	
	  randphase = zeros(size(ptemp,1),duration);
	  
	  for e = 1:size(randphase,1)
	    randphase(e,:) = ptemp(e,tstart(e):tstart(e)+duration-1);
	  end
	  
	  % run the test over time
	  for t = 1:size(randphase,2)
	    [junk,p_ray_boot{i}(loop,f,t)] = rayleigh(double(randphase(:,t)));
	  end
	end
      end
      
      
      if do_disp_anova
	% perform disp_anova test on all pairwise combinations
	count = 0;
	for i = 1:length(events)-1
	  for j = i+1:length(events)
	    count = count+1;
	    
	    % randomize event compairsons for the disp_anova
	    randphase_all = [squeeze(phase{i}(:,f,(1+shift):(duration+shift-1))) ; squeeze(phase{j}(:,f,(1+shift):(duration+shift-1)))];
	    randphase_all = randphase_all(randperm(size(randphase_all,1)),:);
	    randphase1 = randphase_all(1:size(phase{i},1),:);
	    randphase2 = randphase_all(size(phase{i},1)+1:end,:);
	    
	  
	    % loop over time
	    for t = 1:size(randphase1,2)
	      [disp_diff_boot{1,count}(loop,f,t),disp_diff_boot{2,count}(loop,f,t),direction] = disp_anova(double(randphase1(:,t)),double(randphase2(:,t)));
	      %disp_diff_boot{1,count}(loop,f,t) = single(double(disp_diff_boot{1,count}(loop,f,t))*direction);
	      disp_diff_boot{1,count}(loop,f,t) = single( (1 + double(disp_diff_boot{1,count}(loop,f,t))*direction) + ((-direction - 1)/2) );
	    end
	  end
	end
      end
    end
    
    fprintf('\nFreqs(%d): ',length(freqs));
  end
end

fprintf('\n');


