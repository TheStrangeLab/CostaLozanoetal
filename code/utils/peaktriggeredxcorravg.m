function [avg1, avg2, r, l, lags, all_loc] = peaktriggeredxcorravg(sig1,sig2,filt,win,minpeakdistance)
%
%  Event-related peak-triggered averaging (Bragin et al, 1995 J Neurosci)
%
%   INPUTS:
%    sig              - unfiltered signal
%    filt             - bandpassed filtered signal: here we will find the peaks of the wave
%    win              - time window around the peak (in samples)
%    minpeakdistance  - minimum distance between peaks (in samples)
%
%   OUTPUTS:
%    avg    = average of the unfiltered data located around each peak
%    lags   = time samples around each peak
avg1 = [];
avg2 = [];
all_loc={};
r = [];
for k=1:size(sig1,1)
  sig1(k,:) = detrend(sig1(k,:));
  sig2(k,:) = detrend(sig2(k,:));
  [pks, loc] = findpeaks(filt(k,:), 'minpeakdistance', minpeakdistance);
  
  iloc = find(loc > win & loc < length(sig1)-win);
  lags = (-win:win);
  loc = loc(iloc);
  
  [tmpavg1, tmpavg2] = deal(zeros(length(loc), 2*win+1));
  rtmp =[];
  for i=1:length(loc)
    tmpavg1(i,:) = sig1(k,loc(i)-win:loc(i)+win);
    tmpavg2(i,:) = sig2(k,loc(i)-win:loc(i)+win);
    [rtmp(i,:),l]=xcorr(tmpavg1(i,:),tmpavg2(i,:),'biased');
   end
  avg1 = [avg1; tmpavg1];
  avg2 = [avg2; tmpavg2];
  r = [r; rtmp];
  all_loc{k} = loc;
  loc=[];
end

