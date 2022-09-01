function [varargout]=ga_erp(chan,events,subj,erpdurationMS,erpoffsetMS,erpbufferMS,relativeMS,filtfreq,filttype,filtorder,resampledRate,subjfield)
%GA_ERP - calculate a grand average event related potential from EEG data.
%
% Calculate a Grand Average ERP across multiple subjects with 
% 95% confidence intervals for a single channel.  The ERP is
% calculated for each subject and then averaged across subjects.
%
% You can group multiple sessions for a single subject by nesting
% the cell arrays:
%   subj = {'501',{'502','503','504'},'505'};
% would give you three subjects with the second subject having
% three sessions combined when calculating the ERP.
%
% In addition, you can specify different channels to combine across
% subjects, which is helpful for intracranial data where the
% channel numbers do not overlap.  The number of channels must
% match the number of subjects.
%
% FUNCTION:
%   [varargout]=ga_erp(chan,events,subj,erpdurationMS,erpoffsetMS,erpbufferMS,relativeMS,filtfreq,filttype,filtorder)
%
% INPUT ARGS:
%   chan = [44,88];     % electrode #s
%   events = events;    % event data
%   subj = {'BR500','BR501'} % Subject sessions
%   erpdurationMS = 1200;  % erp time length in samples
%   erpoffsetMS = -200;    % offset at which to start the erp in samples
%   erpbufferMS = 1000;    % buffer around each trace (used for
%                       %   filtering, default is 0 for no buffer)
%   relativeMS = [-200 0]  % relative baseline range, 0 for no baseline subtraction
%   filtfreq = [40];    % Filter freq (depends on type, see buttfilt)
%                       %   default is []
%   filttype = 'low';   % Filter type (see buttfilt)
%   filtorder = 1;      % Filter order (see buttfilt)
%   resampledRate = 256; % Resample all the data to this rate
%   subjfield = 'subject' % Name of the event field containing the subj identifier
%
% OUTPUT ARGS:
%   across_sub_ave- ERP averaging across each subject's ERP
%   across_sub_err- 95% Confidence Interval (CI)
%   each_sub_ave- ERP for each subject
%   each_sub_err- 95% CI
%   one_sub_ave- ERP treating all the data as one subject
%   one_sub_err- 95% CI
%   trange- Time range
%

if nargin < 12
  subjfield = 'subject';
  if nargin < 11
    resampledRate = [];
    if nargin < 10
      filtorder = 1;  % first order filter
      if nargin < 9
	filttype = 'low';  % low pass filter
	if nargin < 8
	  filtfreq = [];  % don't filter
	  if nargin < 7
	    relative = 0; % don't subtrace a relative baseline
	    if nargin < 6 
	      eegbuffer = 0; % don't add a buffer
	    end; 
	  end
	end
      end
    end
  end
end

% get initial data info
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(1));

if isempty(resampledRate)
  resampledRate = samplerate;
end

% convert the MS to samples
erpduration = round(erpdurationMS*resampledRate/1000);
erpoffset = round(erpoffsetMS*resampledRate/1000);
erpbuffer = round(erpbufferMS*resampledRate/1000);
relative = round(relativeMS*resampledRate/1000);

if(length(relative) == 2) 
    % adjust the relrange for the offset
    relative = relative - erpoffset+1;
end    

% calculate the ERP for each subject
each_sub_ave = zeros(length(subj),erpduration);
each_sub_err = zeros(length(subj),erpduration);

% for combined subject
g_sum = zeros(1,erpduration);
g_sumsq = zeros(1,erpduration);
g_count = zeros(1,1);


for s = 1:length(subj)
  % status
  fprintf('Subj(%d): ',s);
  
  % generate the string to filter events
  evtstr = '';
  % see if has multiple sessions
  if iscell(subj{s})
    % subj has multiple sessions, so get those events
    for i = 1:length(subj{s})
      if i > 1
	% add the or
	evtstr = [evtstr ' | '];
      end
      
      % see if string
      if isstr(subj{s}{i})
	fprintf('%s ',subj{s}{i});
            
	evtstr = [evtstr 'strcmp(' subjfield ',''' subj{s}{i} ''')'];
      else
	fprintf('%d ',subj{s}{i});
	evtstr = [evtstr subjfield ' == ' num2str(subj{s}{i})];
      end
    end
  else
    
    % see if string
    if isstr(subj{s})
      % status
      fprintf('%s ',subj{s});
      
      % is single session
      evtstr = ['strcmp(' subjfield ',''' subj{s} ''')']; 
    else
	fprintf('%d ',subj{s});
	evtstr = [evtstr subjfield ' == ' num2str(subj{s})];
    end
  end
  
  % get the events for that subject
  s_events = filterStruct(events,evtstr);
  
  % see if there is a special channel
  if length(chan) > 1
    cur_chan = chan(s);
  else
    cur_chan = chan;
  end
  
  % get the data
  eeg=gete_ms(cur_chan,s_events,erpdurationMS,erpoffsetMS,erpbufferMS,filtfreq,filttype,filtorder,resampledRate);
  fprintf('\n');
  
  % see if perform baseline subtraction
  if length(relative) == 2
    % calculate the relative
    releeg = mean(eeg(:,relative(1):relative(2)),2);
    
    % subtract the baseline
    for e=1:length(s_events)
      eeg(e,:)=eeg(e,:)-releeg(e);             
    end;
    
  end 
  
  % do the combined
  g_sum = g_sum + sum(eeg);
  g_sumsq = g_sumsq + sum(eeg.^2);
  g_count = g_count + size(eeg,1);
  
  % calculate the subject's erp
  if length(s_events) > 1
    each_sub_ave(s,:) = mean(eeg);
    each_sub_err(s,:) = 1.96*std(eeg)/sqrt(length(s_events)-1);
  else
    each_sub_ave(s,:) = eeg;
    each_sub_err(s,:) = zeros(size(eeg));
  end
end
  
% get the grand average and error
if length(subj) > 1
  across_sub_ave = mean(each_sub_ave);
  across_sub_err = 1.96*std(each_sub_ave)/sqrt(length(subj)-1);
else
  across_sub_ave = each_sub_ave;
  across_sub_err = zeros(size(eeg));
end


% average and error for combined as single subject
if g_count > 0
  one_sub_ave = g_sum / g_count;
  one_sub_err = sqrt((g_sumsq-g_sum/g_count)/(g_count-1))/sqrt(g_count);
else
  one_sub_ave = 0;
  one_sub_err = 0;
end
  
  
% set the output
for i=1:nargout
  switch i
   case 1
    varargout(i) = {across_sub_ave};
   case 2
    varargout(i) = {across_sub_err};
   case 3
    varargout(i) = {each_sub_ave};
   case 4
    varargout(i) = {each_sub_err};
   case 5
    varargout(i) = {one_sub_ave};
   case 6
    varargout(i) = {one_sub_err};
   case 7
    varargout(i) = {erpoffsetMS:1000/resampledRate:erpdurationMS+erpoffsetMS-1};
  end
end


