function [avg, lags, loc] = peaktriggeredavg(sig,filt,win,minpeakdistance)
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

[pks, loc] = findpeaks(filt, 'minpeakdistance', minpeakdistance);

iloc = find(loc > win & loc < length(sig)-win);
lags = (-win:win);
loc = loc(iloc);

avg = zeros(length(loc), 2*win+1);
for i=1:length(loc)
  avg(i,:) = sig(loc(i)-win:loc(i)+win);
end

