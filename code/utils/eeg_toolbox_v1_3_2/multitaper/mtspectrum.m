function [sp,f,lvar,dof]=mtspectrum(X,tbw,Fs,filtFreq)
%MTSPECTRUM - Calculate Multitaper Spectrum
% 
% MTSPECTRUM calculates the direct spectral estimate of
% series using multitaper methods. Jacknife estimates of the
%        coherence and the tolerances are returned
%
% FUNCTION:
%   [sp,f,lvar,dof]=mtspectrum(X,tbw,Fs,filtFreq)
%
% INPUT ARGS:
%   X: the time series (preferably detrended and
%        standardized.
%   tbw: time-bandwidth factor of the tapers
%   Fs: sampling frequency
%   filtFreq: Single frequency to notch filt using linefilt
%
% OUTPUT ARGS:
%   sp: spectral estimate
%   f: frequencies
%   lvar: jacknife variance.
%   dof: degrees of freedom
%
%

if nargin < 4
  filtFreq = [];
end

% set taper width to sampling rate
width = Fs;
K = 2*tbw-1;

% set number of frequency bins and sampling freq. to 
% next power of 2 above the sampling frequency
nf = 2^nextpow2(Fs);
sampling = nf;

% calc the dpss
E = dpss(width,tbw,'calc');

% first extend the data to nearest second
X = extend(X,Fs);

% make sure the data passed in is in the right dimensions
szX = size(X);
if prod(szX) == length(X)
  % is vector
  N = length(X);
  X = X(:);
  ntrial = 1;
else
  N = szX(2);
  ntrial = szX(1);
  X = X';  
end

% see if we are to filter
if ~isempty(filtFreq)
  for i = 1:ntrial
    X(:,i) = linefilt(X(:,i),filtFreq,Fs,tbw,E);
  end
end

% set the taper overlap and number of windows
M = width;
overlap = M/2;
nwin = fix((N - overlap)/(M-overlap));
if(nwin < 1) 
  nwin=1;
end

% set the frequencies calculated
f = [0:nf/2-1]*sampling/nf;
sp = zeros(ntrial,length(f));
lvar = zeros(ntrial,length(f));
dof = zeros(1,ntrial);

% calculate the spectrum
Xall = X;
for i = 1:ntrial
  X = Xall(:,i);
  y1k = zeros(nf,(K)*nwin);
  tdof = K*nwin;
  for wind = 1:nwin
    xseg = X(overlap*(wind-1)+1:(wind-1)*overlap+M);
    xseg = (xseg-mean(xseg))/std(xseg);
    y1k(:,K*(wind-1)+1:K*wind)=fft(xseg(:,ones(K,1)).*E(:,1:K),nf);
  end;

  %full = y1k(1:nf/2,:).*conj(y1k(1:nf/2,:));

  % compute jacknife estimates and tolerances
  spdel1 = zeros(nf/2,tdof);
  for k =1:tdof
    indices = setdiff([1:tdof],k);
    sp1=sum(y1k(:,indices).*conj(y1k(:,indices)),2);
    spdel1(:,k) = sp1(1:nf/2);
  end;
  
  % now calculate the estimate (Thomson and Chave, 1991 (eqn. 2.3))
  
  xfspdel1 = log(spdel1/(tdof-1)); % do a variance stabilizing log transform
  sp(i,:) = mean(xfspdel1,2)';

  keyboard

  % now estimate the jacknife variance (Thomson and Chave, 1991 (eqn
  % 2.5)) 
  lvar(i,:) = sqrt(tdof-1)*std(xfspdel1,1,2)';

  % get the degrees of freedom
  for k=1:K
    a = sum(E(1:M/2,k).*E(M/2+1:M,k));
    dof(i) = 2*nwin/(1+2*(1-1/nwin)*a^2) + dof(i);
  end;
end

