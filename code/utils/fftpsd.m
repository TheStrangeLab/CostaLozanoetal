
function [psd_avg, faxis]=fftpsd(s, dt, window, pad, normalize)
%
% FFTPSD computes the average power spectrum
% 
% INPUT
%  	s           - signal of interest (1 x N trials)
%  	dt          - sampling interval.
%   window      - number of samples to epoch the data
%   pad         - number of samples to be padded with zeros
%   normalize   - different ways to normalize the single trials before the averaging
%
% OUTPUT
%   psd_avg     - PSD of the signal of interest (1 x N trials)
%   faxis       - frequency axis
%
% 15-Jan-2021 18:07:41, Diego Lozano-Soldevilla CTB-UPM
%

nrpt = floor(size(s,2)/window);
if ~isempty(pad)
  N = pad;
else
  N = window;
end
if isempty(normalize)
  normalize='sum';
end

psd_avg=0;
for k=1:nrpt
  dat = s((k-1)*window+1:k*window);
  power = abs(fft(dat.*hann(length(dat))'-mean(dat.*hann(length(dat))'),pad)).^2;
  power = power(1:N/2+1);
  switch normalize
    case 'sum'
      psd_avg = psd_avg + power;
    case 'sum01'
      psd_avg = psd_avg + (power/max(power));
    case 'db'
      psd_avg = psd_avg + 10.0*log10(power);
    case 'db01'
      psd_avg = psd_avg + 10.0*log10(power / max(power));
  end
end
df = 1.0/(N*dt);
faxis = (0:N/2)*df;
