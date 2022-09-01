function [MI,binAmp] = ModIndex_v2_ft_subsamp(s1,s2,srate,flow,flow_bw,fhigh,fhigh_bw,position,nrpt,nresamp)
% Programmed by Adriano Tort, CBD, BU, 2008 (ModIndex_v1.m). ModIndex_v2_ft
% measures phase-amplitude cross-frequency coupling measure controlling the
% number of trials (See Tort et al J Neurophysiol 2010 and the Supp Info of
% Tort et al PNAS 2008). I only change the bandpass filter function to use
% the Fieldtrip one (DLS)
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
%   nrpt        - number of trials to subsample s1 and s2
%   nresample   - number of time to resample the nrpt trials in s1 and s2
%
% OUTPUT:
%   MI          - modulation index
%   binAmp     - amplitude distribution per phase bin (non-normalized)
%
% 15-Jul-2020 17:11:25 Diego Lozano-Soldevilla
%

nbin=length(position);
winsize = 2*pi/nbin;

if size(s1,1)<nrpt
  error('number of trials is lower than the subsample fold');
end
if size(s1,1)==nrpt
  nresamp=1;
end

mi = zeros(size(flow,2),size(fhigh,2),nrpt);
MI = zeros(size(flow,2),size(fhigh,2),nrpt);
binamp = nan(size(flow,2),size(fhigh,2),nbin,nrpt);
binAmp = nan(size(flow,2),size(fhigh,2),nbin,nrpt);
for n=1:nresamp
  beg=tic;
  rpt = randperm(size(s1,1),nrpt);
  s1_r = s1(rpt,:);
  s2_r = s2(rpt,:);
  for k = 1:size(s1_r,1)
    for x = 1:size(flow,2)
      Fbp_x = [flow(x)-flow_bw/2 flow(x)+flow_bw/2];
      N_x = 3*fix(srate/Fbp_x(1));
      filt_x = ft_preproc_bandpassfilter(zscore(s1_r(k,:)),srate,Fbp_x,N_x,'fir','twopass');
      phs = ft_preproc_hilbert(filt_x,'angle');
      filt_x =[];
      
      for y = 1:size(fhigh,2)
        Fbp_y = [fhigh(y)-fhigh_bw/2 fhigh(y)+fhigh_bw/2];
        N_y = 3*fix(srate/Fbp_y(1));
        filt_y = ft_preproc_bandpassfilter(zscore(s2_r(k,:)),srate,Fbp_y,N_y,'fir','twopass');
        amp = ft_preproc_hilbert(filt_y,'abs');
        filt_y =[];
        
        MeanAmp=nan(1,nbin);
        for j=1:nbin
          phs_mask = find(phs <  position(j)+winsize & phs >=  position(j));
          MeanAmp(1,j)=mean(amp(1,phs_mask));
        end
        phs_mask=[];
        
        mi(x,y,k)=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);
        binamp(x,y,:,k)=MeanAmp;
      end
    end
  end
  MI = nansum(cat(4,MI,mi),4);
  binAmp = nansum(cat(5,binAmp,binamp),5);
  tend=toc(beg);
  disp(['subsample ' num2str(n) ' took ' num2str(tend) 's; ' num2str(nresamp-n) ' to finish']);
end
MI = MI/nresamp;
binAmp = binAmp/nresamp;
