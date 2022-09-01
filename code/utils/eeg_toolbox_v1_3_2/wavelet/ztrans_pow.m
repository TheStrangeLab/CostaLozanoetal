function [zpow,zmean,zstd,kInd,phase] = ztrans_pow(chan,event_info,baseline_info,varargin)
%ZTRANS_POW - Calculate Z-transformed power.
%
% This function will calculate Z-transformed power using the mean and
% std of the power of all events in event_info.  If the baseline_info
% is supplied, it is used to calculate the mean and std, instead.  
%
% The std is calculated across the means of the events (i.e., for a
% particular frequency it calculates the mean of each event and then
% calculates the mean and std across the events).  Thus, it is
% important that you provide enough events to calculate a std.
%
% If you pass in the zmean and zstd as members of the event_info or
% baseline_info structures, it will use those instead of calculating
% them.
%
% This function uses getphasepow to calculate power and converts
% all the power values to singles to save space.  You can pass
% optional args to getphasepow using the varargin at the end.
%
% FUNCTION:
%   [zpow,zmean,zstd,kInd,phase] = ztrans_pow(chan,event_info,baseline_info,varargin)
%
% INPUT ARGS:
%   chan = 42;   % Channel to process
%   event_info = struct('events',{ev},...
%                       'DurationMS',{2000},...
%                       'OffsetMS',{0},...
%                       'BufferMS',{1000});
%   baseline_info = []; % Optional struct defining a baseline
%   varargin -> See getphasepow for optional args.
%
% OUTPUT ARGS:
%   zpow(Events,Freqs,Time) - Z-Transformed Power
%   zmean
%   zstd
%   kInd - Logical indexes of events not thrown out due to kurtosis
%   phase
%
%
% Changes:
%   1/6/06 - PBS - Now returns phase, too, if you want it.
%   07/09/05 MvV - added functionality to see which events were
%      thrown out by kthresh 

if ~exist('baseline_info','var')
  baseline_info = [];
end


% get the event power
events = event_info.events;
DurationMS = event_info.DurationMS;
OffsetMS = event_info.OffsetMS;
BufferMS = event_info.BufferMS;
[phase,zpow,kInd] = getphasepow(chan,events,DurationMS,OffsetMS,BufferMS,'usesingles',varargin{:});
zpow(zpow<=0) = eps;
zpow = log10(zpow);

% get the baseline mean and std for each freq
if ~isempty(baseline_info)
  % get baseline
  % see if passed in already
  if isfield(baseline_info,'zmean') & isfield(baseline_info,'zstd')
    zmean = baseline_info.zmean;
    zstd = baseline_info.zstd;
  else
    % must calculate anew
    events = baseline_info.events;
    DurationMS = baseline_info.DurationMS;
    OffsetMS = baseline_info.OffsetMS;
    BufferMS = baseline_info.BufferMS;
    basepow = getphasepow(chan,events,DurationMS,OffsetMS,BufferMS,'powonly','usesingles',varargin{:});
    basepow(basepow<=0) = eps;
    basepow = log10(basepow);
    
    % get mean and std
    
    % this method calculates the std over time for each event and
    % takes the average std across events
    % (lower freqs often get smaller std)
    zmean = mean(mean(basepow,3),1);
    zstd = mean(std(basepow,0,3),1);
    
    % this method would mean over time in each event then calculate
    % std across events (higher freqs get lower std)
    %zmean = mean(mean(basepow,3),1);
    %zstd = std(mean(basepow,3),0,1);

    % this method would concatenate all the events together and
    % calculate std across all time of all events
    % (equivalent std across freqs)
    %zmean = mean(double(reshape(shiftdim(basepow,2),...
    %                            size(basepow,1)*size(basepow,3),...
    %                            size(basepow,2))'),2);
    %  zstd =  std(double(reshape(shiftdim(basepow,2),...
    %		       size(basepow,1)*size(basepow,3),...
    %		       size(basepow,2))'),0,2);
  end
else
  % use events as baseline
  % Check for zmean ands zstd
  if isfield(event_info,'zmean') & isfield(event_info,'zstd')
    % use provided zmean  
    zmean = event_info.zmean;
    zstd = event_info.zstd;
  else
    zmean = mean(mean(zpow,3),1);
    zstd = mean(std(zpow,0,3),1);
    %zmean = mean(double(reshape(shiftdim(zpow,2),...
    %                            size(zpow,1)*size(zpow,3),...
    %                            size(zpow,2))'),2);
    %  zstd =  std(double(reshape(shiftdim(zpow,2),...
    %		       size(zpow,1)*size(zpow,3),...
    %		       size(zpow,2))'),0,2);
  end
end

% Correct the pow and return the zpow
for f = 1:size(zpow,2)
  zpow(:,f,:) = single((double(zpow(:,f,:)) - zmean(f)) / zstd(f));
end
