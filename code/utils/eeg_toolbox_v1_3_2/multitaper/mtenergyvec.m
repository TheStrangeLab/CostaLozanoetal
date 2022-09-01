function [pow,phase,work] = mtenergyvec(data,freqs,samplerate,bandwidth,windowsize,work)
%MTENERGYVEC - Calculate power and phase with multitapers.
%
% 
% FUNCTION:
%   [pow,phase,work] = mtenergyvec(data,freqs,samplerate,bandwidth,windowsize,work);
%
% INPUT ARGS:
%   data = eeg(1,:);  % A single vector of data to process
%   freqs = [1:60];   % Vector of freqs to process
%   samplerate = 256; % Samplerate of the data
%   bandwidth = 1;    % Frequency bandwith of the tapers
%   windowsize = .3;  % Size of the taper window in seconds
%   work = work;      % Structure returned from previous call.
%                     % This can save processing time for repeated 
%                     % calls on data with the same taper info.
%
% OUTPUT ARGS:
%   pow(freq,time)- Power as a function of time and frequency
%   phase(freq,time)- Phase as a function of time and frequency
%   work- structure of taper and filter info to save for future calls
%
%
%

if ~exist('bandwidth','var')
  bandwidth = 1;
end
if ~exist('windowsize','var')
  windowsize = .3;
end

% set up params
N=fix(samplerate*windowsize); %1 second filter width
K=2*bandwidth-1; % number of tapers
freqs = freqs(:);

% allocate space
pow = zeros(length(freqs),length(data));
phase = zeros(length(freqs),length(data));

% see if must set up tapers and freq filter
if ~exist('work','var') | ...
  isempty(work) | ...
  work.bandwidth~=bandwidth | ...
  work.windowsize~=windowsize | ...
  work.samplerate~=samplerate
  
  % give notification
  fprintf('Setting up new tapers and filter...')
  
  % Set up tapers
  E=dpss(N,bandwidth,K,'calc');
  
  % set up filter for each freq
  pr_op=zeros(length(freqs),N,K);
  x = (-2.*pi.*j.*[1:N]./samplerate);
  shifter = exp(freqs(:,ones(1,length(x))).*x(ones(length(freqs),1),:));
  filt = zeros(length(freqs),2*N);
  
  % apply the shift to each taper for each freq 
  for i=1:K
    pr_op(:,:,i)=shifter.*E(:,i*ones(length(freqs),1))';
  end  
  
  % generate the filter for each freq
  for f = 1:length(freqs)
    a = pr_op(f,:,:);
    x = reshape(a,[size(a,2) size(a,3)]);
    pr = x*x';
    %pr=pr_op(f,:,:)*pr_op(f,:,:)';
    for t=0:N-1 
      filt(f,t+1:t+N)=filt(f,t+1:N+t)+pr(:,N-t)'; 
    end
    filt(f,:)=filt(f,:)*2./N;
    filt(f,:) = filt(f,:)/norm(filt(f,:));
  end

  % get range
  trange = ceil(size(filt,2)/2):(size(filt,2)+length(data)-1)-floor(size(filt,2)/2);
  
  % save the vals to work
  work.filt = filt;
  work.trange = trange;
  work.bandwidth = bandwidth;
  work.windowsize = windowsize;
  work.samplerate = samplerate;
  
  fprintf('Done\n');
else
  % get vars from work
  filt = work.filt;
  trange = work.trange;
end
  
% loop over freqs and calc pow and phase
for f = 1:length(freqs)
  % prepare the filter
  %shifter=exp(-2.*pi.*j.*freqs(f).*[1:N]./samplerate);
  %filt=zeros(1,2.*N);
  %for i=1:K
  %  pr_op(:,i)=shifter'.*E(:,i);
  %end  
  %pr=pr_op*pr_op';
  %for t=0:N-1 
  %  filt(t+1:t+N)=filt(t+1:N+t)+pr(:,N-t)'; 
  %end
  %filt=filt*2./N;

  % convolve the data with the filter
  %y = conv(data,filt/norm(filt));
  y = conv(data,filt(f,:));

  % extract pow and phase
  pow(f,:) = abs(y(trange)).^2;
  
  
  % normalizes phase estimates to length one
  l = find(abs(y) == 0); 
  y(l) = 1;
  y = y./abs(y);
  y(l) = 0;
   
  phase(f,:) = angle( y(trange) );
end

