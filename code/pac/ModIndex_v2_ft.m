function [MI,binAmp] = ModIndex_v2_ft(s1,s2,srate,flow,flow_bw,fhigh,fhigh_bw,position)
% Programmed by Adriano Tort, CBD, BU, 2008 (ModIndex_v1.m). ModIndex_v2_ft 
% measures phase-amplitude cross-frequency coupling measure (See Tort et al 
% J Neurophysiol 2010 and the Supp Info of Tort et al PNAS 2008). I only
% change the bandpass filter function to use the Fieldtrip one (DLS)
%
% INPUT:
%
%  	s1          - signal of interest from channel 1 (1 x N trials)
%  	s2          - signal of interest from channel 2 (1 x N trials)
%   srate       - sampling rate
%   flow        - low frequency axis (1xN vector)
%   flow_bw     - low frequency bandwidth such that [flow(x)-flow_bw/2 flow(x)+flow_bw/2];  (number)
%   flow        - high frequency axis (1xN vector)
%   fhigh_bw    - high frequency bandwidth such that [fhigh(y)-fhigh_bw/2 fhigh(y)+fhigh_bw/2]; (number)
%   position    - low frequency phase binning (number)
%
% OUTPUT:
%   MI          - modulation index
%   binAmp     - amplitude distribution per phase bin (non-normalized)
%
% 21-Sep-2019 10:53:35 Diego Lozano-Soldevilla
%

nbin=length(position);
winsize = 2*pi/nbin;

MI = nan(size(flow,2),size(fhigh,2),size(s1,1));
binAmp = nan(size(flow,2),size(fhigh,2),nbin,size(s1,1));
for k = 1:size(s1,1)
  for x = 1:size(flow,2)
    Fbp_x = [flow(x)-flow_bw/2 flow(x)+flow_bw/2];
    N_x = 3*fix(srate/Fbp_x(1));
    filt_x = ft_preproc_bandpassfilter(zscore(s1(k,:)),srate,Fbp_x,N_x,'fir','twopass');
    phs = ft_preproc_hilbert(filt_x,'angle');
    filt_x =[];
    
    for y = 1:size(fhigh,2)
      Fbp_y = [fhigh(y)-fhigh_bw/2 fhigh(y)+fhigh_bw/2];
      N_y = 3*fix(srate/Fbp_y(1));
      filt_y = ft_preproc_bandpassfilter(zscore(s2(k,:)),srate,Fbp_y,N_y,'fir','twopass');
      amp = ft_preproc_hilbert(filt_y,'abs');
      filt_y =[];
      
      MeanAmp=nan(1,nbin);
      for j=1:nbin
        phs_mask = find(phs <  position(j)+winsize & phs >=  position(j));
        MeanAmp(1,j)=mean(amp(1,phs_mask));
      end
      phs_mask=[];

      MI(x,y,k)=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);
      binAmp(x,y,:,k)=MeanAmp;
    end
  end
end