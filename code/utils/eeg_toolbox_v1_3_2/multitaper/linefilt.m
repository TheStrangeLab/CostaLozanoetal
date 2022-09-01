function yfilt = linefilt(dat,f0,fs,NW,E)
% LINEFILT - Multitaper notch filter.
%
% Filters data set dat removing the line frequency f0
% with the sampling rate of fs using DPSS data tapers
% with time-bandwidth product of NW.  You can either pass in the
% tapers to use to E, or if you leave E as [], this
% function will create the tapers using dpss.
%
% FUNCTION:
%   yfilt = linefilt(dat,f0,fs,NW,E)
%
%

if nargin < 5
  E = [];
end

T = length(dat);
nseg = fix(T/(fs/2));

onesecond=fs;
overlap=onesecond/4;
nseg=(T-onesecond)/overlap+1;

if(fix(nseg)~=nseg) nseg=round(nseg)-1; end;

% get data tapers
K = (2*NW)-1;
if isempty(E)
  E = dpss(onesecond, NW);
end
Wk0 = sum(E(:, 1:K));

fftpad = max(2^nextpow2(onesecond*2),onesecond);
%fftpad=onesecond;

%loop over data segments
for seg = 1:nseg,
  istart=(seg-1)*overlap+1;
  iend=istart+onesecond-1;
  tmp = dat(istart:iend);
  X = fft( E(:,1:K) .* tmp(:, ones(K,1)), fftpad);
  mu(:,seg) = X(1:fftpad/2, :) * Wk0' / sum(Wk0.^2);

end; 

[fstat, muhat, f, trueid] = ftest(dat(1:min(T,512)), NW, K, fftpad, fs, f0);
critval = fix(finv(1-2/fftpad,2,2*K));
trueid=find(fstat>critval);
trueid=trueid(find(trueid>1));
f(trueid);
fstat(trueid);
[fdum,indx]=sort(fstat(trueid));
trueid = trueid(indx);
if(length(indx)>4)
trueid = trueid(end-3:end);
end
%trueid = round(f0*fftpad/fs)+1;
%if(mod(f0,2)) trueid=trueid-1;end;
%phas(1) = 0;
%phas(2:seg) = ph(mu(trueid,1:nseg-1));

phas = zeros(length(trueid),nseg);
phas(:,1:nseg) = ph(mu(trueid,1:nseg));
inum = size(phas,1);

if(nseg > 1)
yfilt=dat;
for ifreq=1:inum
freq = round(f(trueid(ifreq)));
for seg = 1:nseg-1,
istart=(seg-1)*overlap+1;
iend=istart+overlap-1;
yfilt(istart:iend) = yfilt(istart:iend) ...
       - 2*abs(mu(trueid(ifreq),seg))*cos(2*pi*freq*(0:overlap-1)'/fs+phas(ifreq,seg));

end;
end;

tmp=yfilt(iend+1:T);
if(length(tmp)>NW)
  istart=(seg-1)*overlap+1;
  E=dpss(length(tmp),NW,K);
  X=fft(E(:,1:K) .* tmp(:, ones(K,1)), fftpad);
  muend=X(1:round(fftpad/2),:)*Wk0'/sum(Wk0.^2);
  phasend=ph(muend(trueid));
  for ifreq = 1:inum
    freq = round(f(trueid(ifreq)));	
    yfilt(iend+1:T)=yfilt(iend+1:T)-2*abs(muend(trueid(ifreq)))*cos(2*pi*freq*(0:length...
		(tmp)-1)'/fs+phasend(ifreq));
  end;
end;

else

  yfilt=dat;
  X=fft(E(:,1:K) .* tmp(:, ones(K,1)), fftpad);
  muend=X(1:round(fftpad/2),:)*Wk0'/sum(Wk0.^2);
  phasend=ph(muend(trueid));
  for ifreq = 1:inum
    yfilt=yfilt-2*abs(muend(trueid(ifreq)))*cos(2*pi*f(trueid(ifreq))*(0:length...
                (tmp)-1)'/fs+phasend(ifreq));
  end;
end;

yfilt=yfilt(:);
