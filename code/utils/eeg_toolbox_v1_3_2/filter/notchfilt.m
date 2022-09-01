function yfilt = notchfilt(dat,f0,fs,tapers)
%NOTCHFILT - Notch filter with sleppian sequences
%
% yfilt = notchfilt(dat,f0,fs,tapers)
% filters data set dat removing the line frequency f0
% with the sampling rate of fs using DPSS data tapers
% with time-bandwidth product of NW and order K
% tapers are actual dpss tapers or an array
% of the form [NW,K,onesecond]

% keep track of its original size
origLen = length(dat);
if(rem(length(dat), fs) > 0)
  newLen = fs * (floor(length(dat)/fs)+1);
  dat((origLen+1):newLen) = 0;
end;
dat = dat(:);

szT=size(tapers);
if(szT(1)>1)
     E=tapers;
     NW=3;
     K=size(tapers,2);
     onesecond=size(tapers,1);
else
NW=tapers(1);
K=tapers(2);
onesecond=tapers(3);
end;

T = length(dat);
nseg = fix(T/(fs/2));

%onesecond=fs/2;
overlap=onesecond/4;
nseg=(T-onesecond)/overlap+1;
if(fix(nseg)~=nseg) nseg=round(nseg)-1; end;
% get data tapers

Wk0 = sum(E(:, 1:K));

%fftpad = max(2^nextpow2(onesecond),onesecond);
%fftpad=onesecond;
fftpad=4096;

%loop over data segments

for seg = 1:nseg,
  istart=(seg-1)*overlap+1;
  iend=istart+onesecond-1;
  tmp = dat(istart:iend);
  X = fft( E(:,1:K) .* tmp(:, ones(K,1)), fftpad);
  mu(:,seg) = X(1:fftpad/2, :) * Wk0' / sum(Wk0.^2);

end; 


trueid = round(f0*fftpad/fs)+1;
if(mod(f0,2)) trueid=trueid-1;end;
%phas(1:nseg) = ph(mu(trueid,1:nseg));
phas = ph(mu(trueid,1:nseg));


%added code
iend=0;
%end added code

for seg = 1:nseg-1,
  istart=(seg-1)*overlap+1;
  iend=istart+overlap-1;
  sumfilt = zeros(overlap,1);
  for i = 1:length(f0)
    sumfilt =  sumfilt +2*abs(mu(trueid(i),seg))*cos(2*pi*f0(i)*(0:overlap-1)'/fs+phas(i,seg));
  end
  %size(2*abs(mu(trueid,seg)))
  %size(cos((2*pi*f0'*(0:overlap-1))'/fs+(ones(overlap,1)*phas(:,seg)'))')
  %sumfilt = sum(2*abs(mu(trueid,seg))*cos((2*pi*f0'*(0:overlap-1))'/fs+(ones(overlap,1)*phas(:,seg)'))');
  %yfilt(istart:iend) = dat(istart:iend) ...
  %    - 2*abs(mu(trueid,seg))*cos(2*pi*f0*(0:overlap-1)'/fs+phas(seg));
  yfilt(istart:iend) = dat(istart:iend) - sumfilt;

end;

tmp=dat(iend+1:T);
istart=(seg-1)*overlap+1;

X=fft(E(:,1:K) .* tmp(:, ones(K,1)), fftpad);
muend=X(1:round(fftpad/2),:)*Wk0'/sum(Wk0.^2);
phasend=ph(muend(trueid));

sumfilt = zeros(length(tmp),1);

for i = 1:length(f0)
  sumfilt = sumfilt + 2*abs(muend(trueid(i)))*cos(2*pi*f0(i)*(0:length(tmp)-1)'/fs+phasend(i));
end
%sumfilt = sum(2*abs(muend(trueid))*cos(2*pi*f0*(0:length(tmp)-1)'/fs+phasend));
%yfilt(iend+1:T)=dat(iend+1:T)-2*abs(muend(trueid))*cos(2*pi*f0*(0:length...
%		(tmp)-1)'/fs+phasend);
yfilt(iend+1:T)=dat(iend+1:T)-sumfilt;


% restore its original size
yfilt=yfilt(1:origLen);