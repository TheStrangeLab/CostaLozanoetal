function [ifMatrix]=calcIFMatrix(f,amp)
%function ifMatrix=blah(modes,sr)
%This function computes a full (huge) matrix of the entire Hilbert spectra.
%
%INPUT:
% f - freq of each mode (each row is one mode)
% amp - amplitude of each mode. (each row is one mode)
%OUTPUT:
% ifMatrix(freq,timepoint,modeNumber) - boolean matrix indicating whether each
%  mode has an IMF at a particular frequency/timepoint.


freqBins=eeganalparams('freqsEMD');
if exist('amp','var')
  doAmp=true;
  ifMatrix=zeros(length(freqBins),size(f,2),size(f,1),'single');
else
  doAmp=false;
  ifMatrix=false(length(freqBins),size(f,2),size(f,1));
end


for modeNum=1:size(f,1)
  [n,freqNumbers]=histc(f(modeNum,:),freqBins);
  goodTimepoints=find(freqNumbers~=0);
  ind=sub2ind(size(ifMatrix),freqNumbers(goodTimepoints),goodTimepoints,repmat(modeNum,1,length(goodTimepoints)));
  if doAmp
    ifMatrix(ind)= amp(modeNum,goodTimepoints);
  else
    ifMatrix(ind)=true;
  end
end
