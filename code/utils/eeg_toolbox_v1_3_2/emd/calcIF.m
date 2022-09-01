function [f,amp,phase]=calcIF(modes,sr)
%function [f,amp,phase]=calcIF(modes,sr)
%INPUT:
% modes is a matrix where each row is one IMF 
% sr is the sampling rate
%OUTPUT:
% f is the instantaneous frequency for each mode (same size as modes)
% amp is the amplitude for each mode (same size as modes)
% phase is the phase of each mode (same size as modes)


amp=zeros(size(modes));
phase=zeros(size(modes));
t = 2:size(modes,2)-1; 
f = zeros(size(modes));
for m=1:size(modes,1)
  h=hilbert(modes(m,:));
  amp(m,:)=abs(h);
  phase(m,:)=angle(h);
    f(m,:) = [nan 0.5*(angle(-h(t+1).*conj(h(t-1)))+pi)/(2*pi) * sr nan];
end

%f=diff(unwrap(phase(:,[1 1:end]),[],2),1,2)/(2*pi)*sr;

