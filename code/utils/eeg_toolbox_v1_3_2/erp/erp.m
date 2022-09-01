function [ave,err,eeg,kInd]=erp(chan,events,erpduration,erpoffset,erpbuffer,relative,filtfreq,filttype,filtorder,kthresh,resampledrate)
%ERP - calculate an event related potential from EEG data.
%
% Calculate an ERP with 95% confidence intervals for a single channel.
%
% FUNCTION:
%   [ave,err,eeg]=erp(chan,events,erpduration,erpoffset,erpbuffer,relative,filtfreq,filttype,filtorder,kthresh,resampledrate)
%
% INPUT ARGS:
%   chan = 44;          % electrode #
%   events = events;    % event data
%   erpduration = 512;  % erp time length in samples
%   erpoffset = -64;    % offset at which to start the erp in samples
%   erpbuffer = 256;    % buffer around each trace (used for
%                       %   filtering, default is 0 for no buffer)
%   relative = [-64 0]  % relative baseline range, 0 for no
%   baseline subtraction
%   filtfreq = [30];    % Filter freq (depends on type, see buttfilt)
%                       %   default is []
%   filttype = 'low';   % Filter type (see buttfilt)
%   filtorder = 1;      % Filter order (see buttfilt)
%   kthresh = 5         % kurtosis threshold
%   resampledrate = 200 % sample rate after resampling
%
% OUTPUT ARGS:
%   ave- The average ERP for the events
%   err- The 95% confidence intervals for the ERP
%   eeg- The eeg data from each trace to make your own average
%   kInd- if kthresh was used, this indexes the events that were
%   thrown out
%

if nargin < 11
  resampledrate = [];
if nargin < 10
  kthresh = [];
if nargin < 9
  filtorder = 1;  % first order filter
  if nargin < 8
    filttype = 'low';  % low pass filter
    if nargin < 7
      filtfreq = [];  % don't filter
      if nargin < 6
	relative = 0; % don't subtrace a relative baseline
	if nargin < 5 
	  eegbuffer = 0; % don't add a buffer
	end; 
      end
    end
  end
end
end
end

if(length(relative) == 2) 
    % adjust the relrange for the offset
    relative = relative - erpoffset+1;
end    

% get the data
eeg=gete(chan,events,erpduration,erpoffset,erpbuffer,filtfreq,filttype,filtorder,resampledrate);



% see if throw out events with weird kurtosis
if ~isempty(kthresh)
  startsize = size(eeg,1);
  k = kurtosis(eeg,1,2); % do kurtosis in the time dimension (i.e.,
                         % dimension 2)
  goodInd = find(k<=kthresh);
  kInd = setdiff(1:size(eeg,1),goodInd);
  eeg = eeg(goodInd,:);
  sizediff = startsize - size(eeg,1);
  if sizediff > 0
    fprintf('Threw out %d events due to kurtosis...\n',sizediff);
  end
end

% see if perform baseline subtraction
if length(relative) == 2
  % calculate the relative
  releeg = mean(eeg(:,relative(1):relative(2)),2);
  
  % subtract the baseline
  for e=1:size(eeg,1)
    eeg(e,:)=eeg(e,:)-releeg(e);             
  end;
  
end 

if length(events) > 0
  ave = mean(eeg);
  err = 1.96*std(eeg)/sqrt(length(events));
else
  ave = eeg;
  err = zeros(size(eeg));
end


