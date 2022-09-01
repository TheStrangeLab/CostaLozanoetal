function [phase,pow]=multiphasevec2(f,S,Fs,width)
% FUNCTION [phase,pow]=multiphasevec(f,S,Fs,width)
%
% Returns the phase and power as a function of time for a range of
% frequencies (f).
%
% Simply calls phasevec in a loop.
%
% INPUT ARGS:
%   f = [2 4 8];   % Frequencies of interest
%   S = dat;       % Signal to process
%   Fs = 256;      % Sampling frequency
%   width = 6;     % Width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   phase- Phase data [freqs,time]
%   power- Power data [freqs,time]
%
%1/2008: josh modified this function to only compute the fft of the
%original signal once, thus speeding the entire algorithm. now this
%one function is basically a combination of multiphasevec,
%phasevec, and fconv.m


pow = zeros(length(f),length(S));
phase = zeros(length(f),length(S));

dt = 1/Fs;
Ly2=nan; Sfft=nan;

for a=1:length(f)
  st = 1/(2*pi*(f(a)/width)); %that last set of parenthesis is the
                              %most annoying matlab bug ever...
  t=-3.5*st:dt:3.5*st;
  curWave=morlet(f(a),t,width);
  
  Ly=length(S)+length(curWave)-1;
  Ly2=pow2(nextpow2(Ly));
  
  if length(Sfft)~=Ly2 %this line is the trick!!!
    Sfft=fft(S,Ly2); 
  end

  curWaveFFT=fft(curWave,Ly2);
  Y=Sfft.*curWaveFFT; %the heart of fconv
  y=ifft(Y,Ly2);


  y=y(1:1:Ly);

  %all the work is now done. just convert y to a phase and power
  thisPow = abs(y).^2;
  startIndex=ceil(length(t)/2);
  endIndex=length(thisPow)-floor(length(t)/2);
  pow(a,:) = thisPow(startIndex:endIndex);
  
  l = find(abs(y) == 0 );
  y(l) = 1;
  % normalizes phase estimates to length one
  y = y./abs(y);
  y(l) = 0;
  phase(a,:)= angle( y(startIndex:endIndex));

end


