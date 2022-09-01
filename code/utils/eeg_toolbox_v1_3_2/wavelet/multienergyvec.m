function [B,t,f]=multienergyvec(S,f,Fs,width)
% function [B,t,f]=multienergyvec(S,f,Fs,width)
% s : signal
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).

B = zeros(length(f),length(S));

for a=1:length(f)
  B(a,:)=energyvec(f(a),S,Fs,width);
end

t = (1:size(S,2))/Fs;  
