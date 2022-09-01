function [fstat, mu, f, trueid] = ftest(dat, NW, K, fftpad, Fs, f0);
%
% calculate F-statistic, that is, F variance-ratio test.
% F has an F distribution, F(2, 2K-2) degrees of freedom.
% hypothesis test for sinusoidal signals under the assuption of a locally white
% background, which implyies the noise is iid complex Gaussian variables.
% The statistical significant level is chosen to be 1-1/N,
%  N - N Raleigh frequencies
% For details, see Thomson(1987).
%
% Usage:
%   [fstat, mu, f, trueid] = ftest(dat, NW, K, fftpad, Fs, f0);
% Input:
%   dat: time series, a column vector (T x 1)
%   NW:  half-time bandwidth product
%   K:  the numbers of data tapers used, e.g. 2*NW-1
%   fftpad: FFT zero-padding length
%   Fs: sampling rate
%   f0: the freq to be removed
% Output:
%   fstat: F-statistic
%   mu:    estimated amplitude of harmonic signals

%   f:  frequency range
%   trueid: theoretical value of frequence index of line noise
% Note:
%  Utilizing the fact that the taper sum is near zero for antisymmetric tapers
%  should greatly speed the computation up.
%



if nargin < 6
  f0 = 60;
end

if nargin < 5
  Fs = 1;
end

if nargin < 4
  fftpad = 2^nextpow2(length(dat));            
end

dat = dat(:);
T = length(dat);
half_pad = fftpad/2;
f = (0:half_pad-1)*Fs/fftpad;

% trueid = 60*fftpad/Fs;
trueid = f0*fftpad/Fs;
%if trueid ~= round(trueid)
%  error('Theoretical index value of line noise should be an integer !');
%end

E = dpss(T, NW);
Wk0 = sum(E(:, 1:K));

X = fft( E(:,1:K) .* dat(:, ones(1, K)), fftpad);

mu = X(1:half_pad, :) * Wk0' / sum(Wk0.^2);

denom = sum(abs(X(1:half_pad, :) - mu*Wk0).^2, 2);
fstat = (K-1)*abs(mu).^2 * sum(Wk0.^2) ./ denom;

