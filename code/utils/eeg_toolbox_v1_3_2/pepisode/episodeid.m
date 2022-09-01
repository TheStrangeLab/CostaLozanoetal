% [P,holes,detected]=episodeid(timecourse,powerthresh,durthresh,shoulder,plotit)
%
% This script identifies oscillatory episodes where the threshold is given
% leadno can be a vector of leads to analyse.
%
% Parameters:
% timecourse - the wavelet transformed signal to analyse
% powerthresh - wavelet power threshold
% durthresh- minimum time (duration threshold)
% shoulder- optionally, you can include a shoulder so that episodes
%           identification won't be subject to edge artifacts. This value
%           should be specified in _samples_
% [plotit]- optionally set this to 1 to have the script show you how it's
%           working.
%              Red - the wavelet power timecourse you passed in
%                    'timecourse'
%              Green - the 'detected' variable (times where episodes were
%                      found)
%              Cyan - the shoulder (which are excluded from the Pepisode
%                     calculation, but included while detecting episodes)
%
% Returns:
% P- Pepisode (percentage of time occupied by oscillatory episodes)
% holes- start and end times of detected episodes. It has units of
%        seconds.
% detected- a binary vector containing a 1 where an episode was detected
%           and a 0 otherwise

function detected=episodeid(timecourse,powerthresh,durthresh,shoulder,plotit)


TIMECOURSE=timecourse;
if(nargin<5) plotit=0; end
NAN_CODE=-99; samplerate=eegparams('samplerate');
t=(1:length(timecourse))/samplerate;
THRESH=powerthresh*ones(size(timecourse));
timecourse(find(timecourse<THRESH))=0; % zero anything under the threshold, same threshold

detected=zeros(1,length(timecourse)); % to store all episodes

x=(timecourse>0); % a binary vector: pass or not

dx=diff(x); pos=find(dx==1)+1; neg=find(dx==-1)+1; % show the +1 and -1 edges
clear dx;
% now do all the special cases to handle the edges
if(isempty(pos) & isempty(neg))
  if(find(x)>0) H=[1;length(timecourse)]; else H=[]; end % all episode or none
elseif(isempty(pos)) H=[1;neg]; % i.e., starts on an episode, then stops
elseif(isempty(neg)) H=[pos;length(timecourse)]; % starts, then ends on an ep.
else
  if(pos(1)>neg(1))
      pos=[1 pos];
  end; % we start with an episode
  if(neg(length(neg))<pos(length(pos)))
      neg=[neg length(timecourse)];
  end; % we end with an episode
  H=[pos;neg]; % NOTE: by this time, length(pos)==length(neg), necessarily
end; % special-casing, making the H double-vector
clear x pos neg;
if(~isempty(H)) % more than one "hole"
  % find epochs lasting longer than minNcycles*period
  goodep=find((H(2,:)-H(1,:))>=durthresh);
  if(isempty(goodep))
  	H=[];
  else
  	H=H(:,goodep);
	% fprintf(1,'ge ');
  end;
  % this becomes detected vector
  for h=1:size(H,2) detected(H(1,h):H(2,h))=1; end;
end % more than one "hole"

if(shoulder~=0) % handle shoulder
  detected(1:shoulder)=0;
  detected((length(detected)-shoulder):length(detected))=0;
end

% now, consolidate all the intervals

x=(detected>0); % a binary vecor: pass or not
dx=diff(x); pos=find(dx==1)+1; neg=find(dx==-1)+1; % show the +1 and -1 edges
clear dx;
% now do all the special cases to handle the edges
if(isempty(pos) & isempty(neg))
  if(find(x)>0) holes=[1;length(timecourse)]; else holes=[]; end % all episode or none
elseif(isempty(pos)) holes=[1;neg]; % i.e., starts on an episode, then stops
elseif(isempty(neg)) holes=[pos;length(timecourse)]; % starts, then ends on an ep.
else
  if(pos(1)>neg(1))
  	pos=[1 pos]; end; % we start with an episode
  % we end with an episode
  if(neg(length(neg))<pos(length(pos)))
  	neg=[neg length(timecourse)]; end;
  holes=[pos;neg]; % NOTE: by now, length(pos)==length(neg), necessarily
end; % special-casing, making the H double-vector
clear x pos neg;

if(plotit)
  plot((1:length(TIMECOURSE))/samplerate,TIMECOURSE,'r-');
  axis('tight'); replot; ax=axis; hold on;
  miny=min(TIMECOURSE); maxy=max(TIMECOURSE); yval=miny+0.2*(maxy-miny);
  plot(holes/samplerate,yval*ones(size(holes)),'+-g');
  legend('off'); xlabel('Time [s]'); ylabel('Wavelet Power(t)');
  yval=miny+0.1*(maxy-miny);
  plot([1 shoulder]/samplerate,[yval yval],'co-');
  plot((length(TIMECOURSE)-[0 shoulder])/samplerate,[yval yval],'co-');
  hold off;
end % plotit

u=find(detected);
if(isempty(u)), P=0; else, P=length(u)/(length(detected)-2*shoulder); end
