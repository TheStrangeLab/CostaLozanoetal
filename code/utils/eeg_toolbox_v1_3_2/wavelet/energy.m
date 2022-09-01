function [B,t,f] = energy(S,f,Fs,width);
% function [B,t,f] = energy(S,f,Fs,width);
%
% Calculates the average of a time-dependent energy representation of
% multiple trials using a Morlet wavelet method.                            
%
% Input
% -----
% S    : signals = trials x time
% f    : frequencies over which to calculate spectrogram 
% Fs   : sampling frequency
% width: number of cycles in wavelet (> 5 advisable)  
%
% Output
% ------
% t    : time
% f    : frequency
% B    : phase-locking factor = frequency x time
%
% To plot:  imagesc(t,f,B),axis xy; 
%
% See also: PHASEVEC, MORLET, PHASEGRAM  

t = (1:size(S,2))/Fs;  

B = zeros(length(f),size(S,2)); 

for i=1:size(S,1)          
%    fprintf(1,'%d ',i); 
    for j=1:length(f)
        B(j,:) = energyvec(f(j),jdetrend(S(i,:)),Fs,width) + B(j,:);
    end
end
B = B/size(S,1);     
